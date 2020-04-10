#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import json
import os
import sys
import numpy as np
import korali

from epidemics.tools.tools import prepare_folder, save_file, load_file, make_path
from epidemics.tools.compute_credible_intervals import compute_credible_intervals


class EpidemicsBase:

  def __init__( self, fileName=None, defaultProperties={}, **kwargs ):

    self.moduleName = self.__class__.__module__

    defaultProperties = { **defaultProperties,
        'nThreads': 1,
        'silent': False,
        'noSave': False,
        'percentages': [0.5, 0.95, 0.99],
        'dataFolder': './data/'
    }

    if not fileName:

      for x in kwargs:
        if x not in defaultProperties:
          sys.exit(f'\n[epidemics] Unknown input argument: {x}\n')

      for key, value in defaultProperties.items():
        setattr( self, key, kwargs.get(key, value) )

      self.saveInfo ={
        'initials': 'initials.pickle',
        'database': 'data_base.pickle',
        'state': 'state.pickle',
        'korali samples': './_korali_samples/',
        'korali propagation': './_korali_propagation/',
        'inference data': './data_for_inference.pickle',
        'figures': './figures/'
      }

      for x in self.saveInfo.keys():
        self.saveInfo[x] = make_path( *self.save_data_path(), self.saveInfo[x] )

      self.save( fileName=self.saveInfo['initials'] )

      self.data = {}
      self.data['Model'] = {}
      self.data['Propagation'] = {}
      self.data['Validation'] = {}

      self.has_been_called = {
        'sample': False,
        'propagate': False,
        'intervals': False
      }

      self.nParameters = 0
      self.parameters = []

      self.propagatedVariables = {}
      self.standardDeviation = []

      self.credibleIntervals = {}

    else:
      self.load(fileName)

    # these variables cannot be pickled
    self.intervalVariables = []
    self.e = korali.Experiment()




  def save( self, fileName=[] ):
    x = dict(self.__dict__)
    if 'e' in x: del x['e']
    if 'intervalVariables' in x: del x['intervalVariables']
    if fileName==[]: fileName = self.saveInfo['state']
    save_file( x, fileName, 'State', 'pickle')




  def load( self, file ):
    state = load_file( file, 'State', fileType='pickle' )

    for key, value in state.items():
      setattr( self, key, value )




  def get_module_name(self,path):
    return 'epidemics.' + path.split('/epidemics/')[1][:-3].replace( '/','.')

  def get_model_name(self,path):
    return os.path.split(path)[1][:-3]



  def computational_model( s ):
    pass

  def computational_model_propagate( s ):
    pass

  def set_variables_for_interval( s ):
    pass


  def save_data_path( self ):
    return (self.dataFolder,)




  def set_korali_output_files( self, folder):

    self.e['File Output']['Enabled'] = True
    prepare_folder( folder )
    relativeSaveFolder = os.path.relpath(folder, './')
    self.e['File Output']['Path'] = relativeSaveFolder





  def sample( self, nSamples=1000 ):

    self.e = korali.Experiment()

    self.nSamples = nSamples

    self.e['Problem']['Type'] = 'Bayesian/Reference'
    self.e['Problem']['Likelihood Model'] = 'Additive Normal General'
    self.e['Problem']['Reference Data']   = list(self.data['Model']['y-data'])
    self.e['Problem']['Computational Model'] = self.computational_model

    self.e['Solver']['Type'] = 'TMCMC'
    self.e['Solver']['Population Size'] = self.nSamples

    self.set_variables_and_distributions()

    self.set_korali_output_files( self.saveInfo['korali samples'] )

    if(self.silent): e['Console Output']['Verbosity'] = 'Silent'

    self.e['File Output']['Frequency'] = 50
    self.e["Store Sample Information"] = True

    k = korali.Engine()
    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    self.parameters = []
    for j in range(self.nParameters):
      self.parameters.append({})
      self.parameters[j]['Name'] = self.e['Variables'][0]['Name']
      self.parameters[j]['Values'] = np.zeros((self.nSamples,1))
      for k in range(self.nSamples):
        self.parameters[j]['Values'][k] = self.e['Results']['Sample Database'][k][j]

    self.has_been_called['sample'] = True
    self.has_been_called['propagation'] = False
    self.has_been_called['intervals'] = False





  def propagate( self ):

    if not self.has_been_called['sample']:
      print('[Epidemics] Sample before propagation')
      return

    self.e = korali.Experiment()

    self.e['Problem']['Type'] = 'Propagation'
    self.e['Problem']['Execution Model'] = self.computational_model_propagate

    for k in range(self.nParameters):
      self.e['Variables'][k]['Name'] = self.parameters[k]['Name']
      self.e['Variables'][k]['Precomputed Values'] = np.squeeze(self.parameters[k]['Values']).tolist()

    self.e['Solver']['Type'] = 'Executor'

    self.set_korali_output_files(self.saveInfo['korali propagation'])

    if(self.silent):
      self.e['Console Output']['Verbosity'] = 'Silent'

    self.e['Store Sample Information'] = True

    k = korali.Engine()

    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    Ns = self.nSamples
    Nv = self.e['Samples'][0]['Saved Results']['Number of Variables']
    Nt = self.e['Samples'][0]['Saved Results']['Length of Variables']

    varNames = []
    for k in range(Nv):
      varNames.append( self.e['Samples'][0]['Saved Results']['Variables'][k]['Name'] )

    # FIXME: too slow
    self.propagatedVariables = {}
    for i,x in enumerate(varNames):
      self.propagatedVariables[x] = np.zeros((Ns,Nt))
      for k in range(Ns):
        self.propagatedVariables[x][k] = np.asarray( self.e['Samples'][k]['Saved Results']['Variables'][i]['Values'] )

    self.propagatedVariables['Standard Deviation'] = np.zeros((Ns,Nt))
    for k in range(Ns):
      self.propagatedVariables['Standard Deviation'][k] = np.asarray( self.e['Samples'][k]['Saved Results']['Standard Deviation Model'] )

    # TODO clear variable?
    self.e = korali.Experiment()

    self.has_been_called['propagation'] = True
    self.has_been_called['intervals'] = False




  def compute_intervals( self ):

    if not self.has_been_called['propagation']:
      print('[Epidemics] Propagation before intervals')
      return

    self.set_variables_for_interval()

    if( self.silent==False ): print('[Epidemics] Loading files...')

    for x in self.intervalVariables.keys():
      self.intervalVariables[x]['Values'] = self.intervalVariables[x]['Formula'](self.propagatedVariables)

    if( self.silent==False ):
      print('[Epidemics] Finished loading files.')
      print('[Epidemics] Calculate credible intervals...')

    self.credibleIntervals = {}
    for x in self.intervalVariables.keys():
      self.credibleIntervals[x] = compute_credible_intervals( self.intervalVariables[x]['Values'],
                                                              self.propagatedVariables['Standard Deviation'],
                                                              self.percentages );

    if( self.silent==False ): print('[Epidemics] Finished calculating credible intervals...')

    self.has_been_called['intervals'] = True
