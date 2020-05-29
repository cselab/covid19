import numpy as np
import korali

import json
import os
import pickle
import sys
import time
import random
import glob

import matplotlib
import matplotlib.pyplot as plt
plt.ioff()


from epidemics.tools.tools import prepare_folder, make_path, save_file, get_truncated_normal, abort, printlog
from epidemics.tools.compute_credible_intervals import compute_credible_intervals

class EpidemicsBase:

  def __init__( self, **kwargs ):

    self.moduleName = self.__class__.__module__

    self.nThreads    = kwargs.pop('nThreads', 1)
    self.silent      = kwargs.pop('silent', False)
    self.silentPlot  = kwargs.pop('silentPlot', False)
    self.noSave      = kwargs.pop('noSave', False)
    self.dataFolder  = kwargs.pop('dataFolder', './data/')
    self.sampler     = kwargs.pop('sampler','TMCMC')
    self.maxGen      = kwargs.pop('maxGen', 100)
    self.display     = os.environ['HOME']
    self.synthetic   = kwargs.pop('synthetic', False)

    if(self.synthetic):
        self.datafile     = kwargs.pop('dataFile')

    if(not self.display):
        matplotlib.use('Agg')

    if kwargs:
        abort(f"Unknown input arguments: {kwargs}")

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

  def load(fileName="state.pickle"):
    """Load object pickled by `save()`"""
    with open(fileName, 'rb') as f:
      return pickle.load(f)


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
    """
    Returns `tuple`:
    Directories to be joined to generate path for output data.
    """
    return (self.dataFolder,)




  def set_korali_output_files( self, folder, frequency = 1 ):

    self.e['File Output']['Enabled'] = True
    self.e['File Output']['Frequency'] = frequency
    prepare_folder( folder )
    relativeSaveFolder = os.path.relpath(folder, './')
    self.e['File Output']['Path'] = relativeSaveFolder


  def get_uniform_priors( self, *triples ) :

    js = {}
    js['Variables']=[]
    js['Distributions']=[]

    k = 0
    for (varname, lowerbound, upperbound) in triples:
      js['Variables'].append({})
      js['Variables'][k]['Name'] = varname
      js['Variables'][k]['Prior Distribution'] = 'Prior for ' + varname

      js['Distributions'].append({})
      js['Distributions'][k]['Name'] = 'Prior for ' + varname
      js['Distributions'][k]['Type'] = 'Univariate/Uniform'
      js['Distributions'][k]['Minimum'] = lowerbound
      js['Distributions'][k]['Maximum'] = upperbound
      k += 1

    return js


  def sample(self, nSamples=1000, cov=0.4):

    self.e = korali.Experiment()

    self.nSamples = nSamples

    self.e['Problem']['Type'] = 'Bayesian/Reference'
    self.e['Problem']['Likelihood Model'] = self.likelihoodModel
    self.e['Problem']['Reference Data']   = list(map(float, self.data['Model']['y-data']))
    self.e['Problem']['Computational Model'] = self.computational_model

    self.e['Solver']['Type'] = "TMCMC"
    self.e['Solver']['Version'] = self.sampler
    self.e['Solver']['Step Size'] = 0.1
    self.e['Solver']['Population Size'] = self.nSamples
    self.e['Solver']['Target Coefficient Of Variation'] = cov
    self.e['Solver']['Termination Criteria']['Max Generations'] = self.maxGen
    js = self.get_variables_and_distributions()
    self.set_variables_and_distributions(js)

    self.set_korali_output_files( self.saveInfo['korali samples'], 100 )
    self.e['Console Output']['Verbosity'] = 'Detailed'
    if(self.silent): self.e['Console Output']['Verbosity'] = 'Silent'

    k = korali.Engine()
    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    js = {}
    js['Evidence'] = self.e['Solver']['LogEvidence']
    printlog(f"Log Evidence = {js['Evidence']}")
    save_file( js, self.saveInfo['evidence'], 'Log Evidence', fileType='json' )

    printlog('Copy variables from Korali to Epidemics...')
    self.parameters = []
    myDatabase = self.e['Results']['Sample Database']
    for j in range(self.nParameters):
      self.parameters.append({})
      self.parameters[j]['Name'] = self.e['Variables'][0]['Name']
      self.parameters[j]['Values'] = np.asarray( [myDatabase[k][j] for k in range(self.nSamples)] )

    self.has_been_called['sample'] = True
    self.has_been_called['propagate'] = False
    printlog('Done copying variables.')


  def optimize( self, nSamples ):

    self.nSamples = 1

    self.e = korali.Experiment()

    self.e['Problem']['Type'] = 'Bayesian/Reference'
    self.e['Problem']['Likelihood Model'] = self.likelihoodModel
    self.e['Problem']['Reference Data']   = list(map(float, self.data['Model']['y-data']))
    self.e['Problem']['Computational Model'] = self.computational_model

    self.e["Solver"]["Type"] = "CMAES"
    self.e["Solver"]["Population Size"] = nSamples
    self.e["Solver"]["Termination Criteria"]["Max Generations"] = self.maxGen
    self.e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-12

    js = self.get_variables_and_distributions()
    self.set_variables_and_distributions(js)

    self.set_korali_output_files( self.saveInfo['korali samples'], 1 )

    if(self.silent): e['Console Output']['Verbosity'] = 'Silent'

    k = korali.Engine()
    k['Conduit']['Type'] = 'Concurrent'
    k['Conduit']['Concurrent Jobs'] = self.nThreads

    k.run(self.e)

    printlog('Copy variables from Korali to Epidemics...')
    self.parameters = []
    myDatabase = self.e['Results']['Best Sample']['Parameters']
    for j in range(self.nParameters):
      self.parameters.append({})
      self.parameters[j]['Name'] = self.e['Variables'][j]['Name']
      self.parameters[j]['Values'] = np.asarray( [myDatabase[j]] )

    self.has_been_called['optimize'] = True
    self.has_been_called['propagate'] = False
    printlog('Done copying variables.')


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

      printlog('Loading posterior samples')

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

      printlog('Loaded')


  def propagate( self, nPropagate = 1000 ):

    if not self.has_been_called['sample'] and not self.has_been_called['optimize'] :
      abort('[Error] Sample or Optimize before propagation')
      return

    self.e = korali.Experiment()

    self.nPropagate = nPropagate

    self.e['Problem']['Type'] = 'Propagation'
    self.e['Problem']['Execution Model'] = self.computational_model_propagate

    for k in range(self.nParameters):
      self.e['Variables'][k]['Name'] = self.parameters[k]['Name']
      self.e['Variables'][k]['Precomputed Values'] = self.parameters[k]['Values'].tolist()

    self.e['Solver']['Type'] = 'Executor'

    self.set_korali_output_files( self.saveInfo['korali propagation'] )

    if(self.silent): self.e['Console Output']['Verbosity'] = 'Silent'

    self.e['Store Sample Information'] = True

    k = korali.Engine()

    k.run(self.e)

    propagate_idx = random.sample(range(self.nSamples), nPropagate)

    Nv = self.e['Samples'][0]['Saved Results']['Number of Variables']
    Nt = self.e['Samples'][0]['Saved Results']['Length of Variables']
    varNames = []
    for k in range(Nv):
      varNames.append( self.e['Samples'][0]['Saved Results']['Variables'][k]['Name'] )

    printlog('Copy variables from Korali to Epidemics...')
    self.propagatedVariables = {}
    for i,x in enumerate(varNames):
      self.propagatedVariables[x] = np.zeros((nPropagate,Nt))
      for k, idx in enumerate(propagate_idx):
        self.propagatedVariables[x][k] = np.asarray( self.e['Samples'][idx]['Saved Results']['Variables'][i]['Values'] )

    if( self.likelihoodModel=='Normal' or self.likelihoodModel=='Positive Normal' ):
      varName = 'Standard Deviation'
    elif( self.likelihoodModel=='Negative Binomial' ):
      varName = 'Dispersion'
    else:
      abort('Likelihood not found in propagate.')

    self.propagatedVariables[varName] = np.zeros((nPropagate,Nt))
    for k in range(nPropagate):
      self.propagatedVariables[varName][k] = np.asarray( self.e['Samples'][k]['Saved Results'][varName] )

    printlog('Done copying variables.')

    # TODO clear variable?
    self.e = korali.Experiment()

    self.has_been_called['propagate'] = True


  def new_figure(self):
    printlog('New figure...')
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(self.modelDescription + '  (' + self.country + ')')

    if(self.silentPlot): plt.ion()

    return fig

  def compute_intervals(self, varName, ns, percentages, cumsum=False):
    Np = self.propagatedVariables[varName].shape[0]
    Nt = self.propagatedVariables[varName].shape[1]

    samples = np.zeros((Np*ns,Nt))

    printlog(f"Sampling from {self.likelihoodModel} for '{varName}' variable... ")

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
        try:
          x = [ np.random.negative_binomial(r,1-p) for _ in range(ns) ]
        except:
          printlog("Error p: {}".format(p))
        samples[:,k] = np.asarray(x).flatten()

    else:
      raise NotImplementedError(
              "Likelihood not found in compute_intervals.")

    if cumsum:
      samples = np.cumsum(samples, 1)

    elapsed = time.process_time() - start
    printlog(f" elapsed {elapsed:.2f} sec")

    printlog(f"Computing quantiles... ")

    mean   = np.zeros(Nt)
    median = np.zeros(Nt)
    for k in range(Nt):
      median[k] = np.quantile(samples[:,k],0.5)
      mean[k] = np.mean(samples[:,k])

    intervals = []
    for p in np.sort(percentages)[::-1]:
      lo = np.zeros(Nt);
      hi = np.zeros(Nt);
      for k in range(Nt):
        lo[k] = np.quantile(samples[:,k],0.5-p/2)
        hi[k] = np.quantile(samples[:,k],0.5+p/2)
      intervals.append((p, lo, hi))

    t = self.data['Propagation']['x-data']

    return t, mean, median, intervals

  def compute_plot_intervals( self, varName, ns, ax, ylabel, cummulate=-1):

    Np = self.propagatedVariables[varName].shape[0]
    Nt = self.propagatedVariables[varName].shape[1]

    samples = np.zeros((Np*ns,Nt))

    printlog(f"Sampling from {self.likelihoodModel} for '{varName}' variable... ", end='', flush=True)

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
        try:
          x = [ np.random.negative_binomial(r,1-p) for _ in range(ns) ]
        except:
          printlog("Error p: {}".format(p))
        samples[:,k] = np.asarray(x).flatten()

    else:
      abort("Likelihood not found in compute_plot_intervals.")

    if cummulate>0 :
      samples = np.cumsum(samples,axis=cummulate)

    elapsed = time.process_time() - start
    printlog(f" elapsed {elapsed:.2f} sec")

    printlog(f"Computing quantiles... ")

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
    x = range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) )
    ax.set_xticks( x[0:-1:3] )
    ax.grid()
    if( self.logPlot ): ax.set_yscale('log')

    plt.draw()
    plt.pause(0.001)

    return samples

def load_param_samples(datadir):
    """
    Returns parameters and log-likelihood from the last generation of samples.
    datadir: `str`
        Path to directory containing `_korali_samples`
    """
    samplespath = sorted(
        glob.glob(os.path.join(datadir, '_korali_samples', '*.json')))[-1]

    with open(samplespath) as f:
        js = json.load(f)

    names = [v['Name'] for v in js['Variables']]
    names_add = ['logPrior', 'logLikelihood']
    names_all = names + names_add

    samples = np.array(js['Solver']['Sample Database'])
    logprior = np.array(js['Solver']['Sample LogPrior Database'])
    loglike = np.array(js['Solver']['Sample LogLikelihood Database'])

    dtype = np.dtype({
            'names' : names_all,
            'formats' : [np.float] * len(names_all)})
    comb = np.empty(samples.shape[0], dtype=dtype)
    for i,name in enumerate(names):
        comb[name] = samples[:,i]
    comb[names_add[0]] = logprior
    comb[names_add[1]] = loglike
    return comb

