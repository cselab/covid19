#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import numpy as np
import korali

import json
import os
import pickle
import sys

import matplotlib.pyplot as plt
plt.ioff()

from epidemics.tools.tools import prepare_folder, make_path
from epidemics.tools.compute_credible_intervals import compute_credible_intervals


class EpidemicsBase:

  def __init__( self, **kwargs ):

    self.moduleName = self.__class__.__module__

    self.nThreads    = kwargs.pop('nThreads', 1)
    self.silent      = kwargs.pop('silent', False)
    self.noSave      = kwargs.pop('noSave', False)
    self.dataFolder  = kwargs.pop('dataFolder', './data/')

    if kwargs:
        sys.exit(f"\n[Epidemics] Unknown input arguments: {kwargs}\n")

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

    # these variables cannot be pickled
    self.intervalVariables = []
    self.e = None  #korali.Experiment()





  def computational_model( s ):
    pass

  def computational_model_propagate( s ):
    pass

  def set_variables_for_interval( s ):
    pass




  def save( self, fileName=None ):
    """Pickle itself to the given target file."""
    if not fileName:
      fileName = self.saveInfo['state']
    with open(fileName, 'wb') as f:
      pickle.dump(self, f)




  def __getstate__(self):
    """Return the state for pickling."""
    state = self.__dict__.copy()
    if 'e' in state:
      del state['e']
    if 'intervalVariables' in state:
      del state['intervalVariables']
    return state




  def get_module_name(self,path):
    return 'epidemics.' + path.split('/epidemics/')[1][:-3].replace( '/','.')




  def get_model_name(self,path):
    return os.path.split(path)[1][:-3]



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
    self.e['Problem']['Likelihood Model'] = self.likelihoodModel + ' General'
    self.e['Problem']['Reference Data']   = list(map(float, self.data['Model']['y-data']))
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

    # FIXME: too slow
    print('[Epidemics] Copy variables from Korali to Epidemics...')
    self.parameters = []
    for j in range(self.nParameters):
      self.parameters.append({})
      self.parameters[j]['Name'] = self.e['Variables'][0]['Name']
      self.parameters[j]['Values'] = np.zeros((self.nSamples,1))
      for k in range(self.nSamples):
        self.parameters[j]['Values'][k] = self.e['Results']['Sample Database'][k][j]

    self.has_been_called['sample'] = True
    self.has_been_called['propagation'] = False
    print('[Epidemics] Done copying variables.')




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
    print('[Epidemics] Copy variables from Korali to Epidemics...')
    self.propagatedVariables = {}
    for i,x in enumerate(varNames):
      self.propagatedVariables[x] = np.zeros((Ns,Nt))
      for k in range(Ns):
        self.propagatedVariables[x][k] = np.asarray( self.e['Samples'][k]['Saved Results']['Variables'][i]['Values'] )

    if( self.likelihoodModel=='Additive Normal' ):   varName = 'Standard Deviation Model'
    if( self.likelihoodModel=='Negative Binomial' ): varName = 'Dispersion'

    self.propagatedVariables[varName] = np.zeros((Ns,Nt))
    for k in range(Ns):
      self.propagatedVariables[varName][k] = np.asarray( self.e['Samples'][k]['Saved Results'][varName] )

    print('[Epidemics] Done copying variables.')

    # TODO clear variable?
    self.e = korali.Experiment()

    self.has_been_called['propagation'] = True
    self.has_been_called['intervals'] = False




  def new_figure(self):
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(self.modelDescription + '  (' + self.country + ')')
    return fig



  def compute_plot_intervals( self, varName, ns, ax, ylabel, cummulate=-1):

    Ns = self.propagatedVariables[varName].shape[0]
    Nt = self.propagatedVariables[varName].shape[1]

    samples = np.zeros((Ns*ns,Nt))

    if self.likelihoodModel=='Additive Normal':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        r = self.propagatedVariables['Standard Deviation Model'][:,k]
        x = [ np.random.normal(m,r) for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    if self.likelihoodModel=='Negative Binomial':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        r = self.propagatedVariables['Dispersion'][:,k]
        p = p =  m/(m+r)
        x = [ np.random.negative_binomial(r,1-p) for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    if cummulate>0 :
      samples = np.cumsum(samples,axis=cummulate)

    mean   = np.zeros((Nt,1))
    median = np.zeros((Nt,1))
    for k in range(Nt):
      median[k] = np.quantile( samples[:,k],0.5)
      mean[k] = np.mean( samples[:,k] )

    for p in np.sort(self.percentages)[::-1]:
      q1 = np.zeros((Nt,));
      q2 = np.zeros((Nt,));
      for k in range(Nt):
        q1[k] = np.quantile( samples[:,k],0.5-p/2)
        q2[k] = np.quantile( samples[:,k],0.5+p/2)
      ax.fill_between( self.data['Propagation']['x-data'], q1 , q2,  alpha=0.5, label=f' {100*p:.1f}% credible interval' )

    ax.plot( self.data['Propagation']['x-data'], mean, '-', lw=2, label='Mean', color='black')
    ax.plot( self.data['Propagation']['x-data'], median, '--', lw=2, label='Median', color='black')

    ax.legend(loc='upper left')
    ax.set_ylabel( ylabel )
    ax.set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax.grid()
    if( self.logPlot ): ax.set_yscale('log')

    return samples
