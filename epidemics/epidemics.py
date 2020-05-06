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
import time

import matplotlib.pyplot as plt
plt.ioff()

from epidemics.tools.tools import prepare_folder, make_path, save_file, get_truncated_normal
from epidemics.tools.compute_credible_intervals import compute_credible_intervals




class EpidemicsBase:

  def __init__( self, **kwargs ):

    self.moduleName = self.__class__.__module__

    self.nThreads    = kwargs.pop('nThreads', 1)
    self.silent      = kwargs.pop('silent', False)
    self.silentPlot  = kwargs.pop('silentPlot', False)
    self.noSave      = kwargs.pop('noSave', False)
    self.dataFolder  = kwargs.pop('dataFolder', './data/')
    self.backend     = kwargs.pop('backend','numpy')

    if kwargs:
        sys.exit(f"\n[Epidemics] Unknown input arguments: {kwargs}\n")

    self.saveInfo ={
      'state': 'state.pickle',
      'korali samples': './_korali_samples/',
      'korali propagation': './_korali_propagation/',
      'inference data': './data_for_inference.pickle',
      'figures': './figures/',
      'evidence': 'evidence.json'
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
      'optimize': False
    }

    self.nParameters = 0
    self.parameters = []
    self.bestParameters = []

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
    self.e['Problem']['Likelihood Model'] = self.likelihoodModel
    self.e['Problem']['Reference Data']   = list(map(float, self.data['Model']['y-data']))
    self.e['Problem']['Computational Model'] = self.computational_model

    self.e['Solver']['Type'] = 'TMCMC'
    self.e['Solver']['Population Size'] = self.nSamples

    js = self.get_variables_and_distributions()
    self.set_variables_and_distributions(js)

    self.set_korali_output_files( self.saveInfo['korali samples'] )

    if(self.silent): self.e['Console Output']['Verbosity'] = 'Silent'

    self.e['File Output']['Frequency'] = 50
    self.e["Store Sample Information"] = True

    k = korali.Engine()
    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    if self.e['Solver']['Type']=='TMCMC':
      js = {}
      js['Evidence'] = self.e['Solver']['LogEvidence']
      print(f"[Epidemics] Log Evidence = {js['Evidence']}")
      save_file( js, self.saveInfo['evidence'], 'Log Evidence', fileType='json' )

    print('[Epidemics] Copy variables from Korali to Epidemics...')
    self.parameters = []
    myDatabase = self.e['Results']['Sample Database']
    for j in range(self.nParameters):
      self.parameters.append({})
      self.parameters[j]['Name'] = self.e['Variables'][0]['Name']
      self.parameters[j]['Values'] = np.asarray( [myDatabase[k][j] for k in range(self.nSamples)] )


    logP = []
    for k in range(self.nSamples):
      logP.append( self.e['Samples'][k]['logPosterior'] )
    ind = np.argmax(np.asarray(logP))

    self.bestParameters = {}
    self.bestParameters['Reference Evaluations'] = self.e['Samples'][0]['Reference Evaluations']
    self.bestParameters['Sample Index'] = ind
    if( self.likelihoodModel=='Normal' or self.likelihoodModel=='Positive Normal' ):
      self.bestParameters['Standard Deviation'] = self.e['Samples'][0]['Standard Deviation']
    elif( self.likelihoodModel=='Negative Binomial' ):
      self.bestParameters['Dispersion'] = self.e['Samples'][0]['Dispersion']
    else:
      sys.exit('\n[Epidemics] Likelihood not found in sample.\n')

    self.has_been_called['sample'] = True
    self.has_been_called['propagate'] = False
    print('[Epidemics] Done copying variables.')




  def optimize( self, maxGenerations):

    self.nSamples = 1
    self.maxGenerations = maxGenerations

    self.e = korali.Experiment()

    self.e['Problem']['Type'] = 'Bayesian/Reference'
    self.e['Problem']['Likelihood Model'] = self.likelihoodModel
    self.e['Problem']['Reference Data']   = list(map(float, self.data['Model']['y-data']))
    self.e['Problem']['Computational Model'] = self.computational_model

    self.e["Solver"]["Type"] = "CMAES"
    self.e["Solver"]["Population Size"] = 24
    self.e["Solver"]["Termination Criteria"]["Max Generations"] = self.maxGenerations
    self.e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-8

    js = self.get_variables_and_distributions()
    self.set_variables_and_distributions(js)

    self.set_korali_output_files( self.saveInfo['korali samples'] )

    if(self.silent): e['Console Output']['Verbosity'] = 'Silent'

    self.e['File Output']['Frequency'] = 10
    self.e["Store Sample Information"] = True

    k = korali.Engine()
    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    print('[Epidemics] Copy variables from Korali to Epidemics...')
    self.parameters = []
    myDatabase = self.e['Results']['Best Sample']['Parameters']
    for j in range(self.nParameters):
      self.parameters.append({})
      self.parameters[j]['Name'] = self.e['Variables'][j]['Name']
      self.parameters[j]['Values'] = np.asarray( [myDatabase[j]] )

    self.has_been_called['optimize'] = True
    self.has_been_called['propagate'] = False
    print('[Epidemics] Done copying variables.')




  def least_squares( self, maxGenerations):

    self.nSamples = 1
    self.maxGenerations = maxGenerations

    self.e = korali.Experiment()

    self.e['Problem']['Type'] = 'Optimization/Stochastic'
    y = list(map(float, self.data['Model']['y-data']))
    self.e['Problem']['Objective Function'] = lambda p: self.sum_of_squares_errors(p,y)

    self.e["Solver"]["Type"] = "CMAES"
    self.e["Solver"]["Population Size"] = 48
    self.e["Solver"]["Termination Criteria"]["Max Generations"] = self.maxGenerations
    self.e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-8

    js = self.get_variables_and_distributions()
    self.set_variables_and_distributions(js)

    self.set_korali_output_files( self.saveInfo['korali samples'] )

    if(self.silent): self.e['Console Output']['Verbosity'] = 'Silent'
    else: self.e['Console Output']['Verbosity'] = 'Detailed'

    self.e['File Output']['Frequency'] = 10
    self.e["Store Sample Information"] = True

    k = korali.Engine()
    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    return self.e['Results']['Best Sample']['Parameters']




  def sum_of_squares_errors(self,s,y):

    self.computational_model( s )

    dif = np.asarray(y) - np.asarray(s['Reference Evaluations'])

    s['F(x)'] = -np.sum(dif*dif) / len(y)




  def set_variables_and_distributions( self, js ):

    if self.e['Problem']['Type']=='Bayesian/Reference' :
      nP = self.nParameters
      for k in range(nP):
        self.e['Variables'][k]['Name'] = js['Variables'][k]['Name']
        self.e['Variables'][k]['Prior Distribution'] = js['Variables'][k]['Prior Distribution']
        self.e['Distributions'][k]['Name'] = js['Distributions'][k]['Name']
        self.e['Distributions'][k]['Type'] = js['Distributions'][k]['Type']
        self.e['Distributions'][k]['Minimum'] = js['Distributions'][k]['Minimum']
        self.e['Distributions'][k]['Maximum'] = js['Distributions'][k]['Maximum']

    else:
      nP = self.nParameters-1
      for k in range(nP):
        self.e['Variables'][k]['Name'] = js['Variables'][k]['Name']


    if self.e["Solver"]["Type"] == "CMAES":
      for k in range(nP):
        self.e["Variables"][k]["Lower Bound"] = js['Distributions'][k]['Minimum']
        self.e["Variables"][k]["Upper Bound"] = js['Distributions'][k]['Maximum']



  def load_parameters(self,samples_path):

      print('[Epidemics] Loading posterior samples')

      files = list(set([filename for filename in os.listdir(samples_path) if (filename.endswith(".json"))]))
      files.sort()
      filename = files[-1]

      variable_names = []
      with open(samples_path+'/'+files[-1]) as json_file:
        data = json.load(json_file)
        samples = data['Results']['Sample Database']
        variables = data['Variables']

      self.nParameters = len(variables)
      self.nSamples = len(samples)
      self.parameters = []
      for j in range(self.nParameters):
        self.parameters.append({})
        self.parameters[j]['Name'] = variables[j]['Name']
        self.parameters[j]['Values'] = np.asarray( [samples[k][j] for k in range(self.nSamples)] )

      self.has_been_called['sample'] = True

      print('[Epidemics] Loaded')


  def propagate( self ):

    if not self.has_been_called['sample'] and not self.has_been_called['optimize'] :
      sys.exit('[Epidemics] [Error] Sample or Optimize before propagation')
      return

    self.e = korali.Experiment()

    self.e['Problem']['Type'] = 'Propagation'
    self.e['Problem']['Execution Model'] = self.computational_model_propagate

    for k in range(self.nParameters):
      self.e['Variables'][k]['Name'] = self.parameters[k]['Name']
      self.e['Variables'][k]['Precomputed Values'] = self.parameters[k]['Values'].tolist()

    self.e['Solver']['Type'] = 'Executor'

    self.set_korali_output_files(self.saveInfo['korali propagation'])

    if(self.silent): self.e['Console Output']['Verbosity'] = 'Silent'

    self.e['Store Sample Information'] = True

    k = korali.Engine()

    k.run(self.e)

    Ns = self.nSamples
    Nv = self.e['Samples'][0]['Saved Results']['Number of Variables']
    Nt = self.e['Samples'][0]['Saved Results']['Length of Variables']

    varNames = []
    for k in range(Nv):
      varNames.append( self.e['Samples'][0]['Saved Results']['Variables'][k]['Name'] )

    print('[Epidemics] Copy variables from Korali to Epidemics...')
    self.propagatedVariables = {}
    for i,x in enumerate(varNames):
      self.propagatedVariables[x] = np.zeros((Ns,Nt))
      for k in range(Ns):
        self.propagatedVariables[x][k] = np.asarray( self.e['Samples'][k]['Saved Results']['Variables'][i]['Values'] )

    if( self.likelihoodModel=='Normal' or self.likelihoodModel=='Positive Normal' ):
      varName = 'Standard Deviation'
    elif( self.likelihoodModel=='Negative Binomial' ):
      varName = 'Dispersion'
    else:
      sys.exit('\n[Epidemics] Likelihood not found in propagate.\n')

    self.propagatedVariables[varName] = np.zeros((Ns,Nt))
    for k in range(Ns):
      self.propagatedVariables[varName][k] = np.asarray( self.e['Samples'][k]['Saved Results'][varName] )

    print('[Epidemics] Done copying variables.')

    # TODO clear variable?
    self.e = korali.Experiment()

    self.has_been_called['propagate'] = True




  def new_figure(self):
    print('[Epidemics] Compute and Plot credible intervals.')
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(self.modelDescription + '  (' + self.country + ')')

    if(self.silentPlot): plt.ion()

    return fig




  def compute_plot_intervals( self, varName, ns, ax, ylabel, cummulate=-1):

    Ns = self.propagatedVariables[varName].shape[0]
    Nt = self.propagatedVariables[varName].shape[1]

    samples = np.zeros((Ns*ns,Nt))

    print(f"[Epidemics] Sampling from {self.likelihoodModel} for '{varName}' variable... ", end='', flush=True)

    start = time.process_time()

    if self.likelihoodModel=='Normal':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        r = self.propagatedVariables['Standard Deviation'][:,k]
        x = [ np.random.normal(m,r) for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    elif self.likelihoodModel=='Positive Normal':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        s = self.propagatedVariables['Standard Deviation'][:,k]
        t = get_truncated_normal(m,s,0,np.Inf)
        x = [ t.rvs() for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    elif self.likelihoodModel=='Negative Binomial':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        r = self.propagatedVariables['Dispersion'][:,k]
        p =  m/(m+r)
        print(1-p)
        x = [ np.random.negative_binomial(r,1-p) for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    else:
      sys.exit("\n[Epidemics] Likelihood not found in compute_plot_intervals.\n")

    if cummulate>0 :
      samples = np.cumsum(samples,axis=cummulate)

    elapsed = time.process_time() - start
    print(f" elapsed {elapsed:.2f} sec")

    print(f"[Epidemics] Computing quantiles... ")

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
    ax.plot( self.data['Propagation']['x-data'], median, '--', lw=2, label='Median', color='blue')

    ax.legend(loc='upper left')
    ax.set_ylabel( ylabel )
    ax.set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax.grid()
    if( self.logPlot ): ax.set_yscale('log')

    plt.draw()
    plt.pause(0.001)

    return samples
