import numpy as np
from epidemics.tools.tools import prepare_folder, save_file
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.sir_int.tnrm'
    self.modelDescription = 'Fit SIR with interventions on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

    self.process_data()


  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )

 
  def process_data( self ):
    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize
    I0 = y[0]
    S0 = N - I0
    y0 = S0, I0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = np.diff( y[0:] )
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = np.diff( y[0:-self.nValidation] )
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = np.diff( y[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )



  def get_variables_and_distributions( self ):

    p = ['beta','gamma', 'tact', 'dtact', 'kbeta', 'Sigma']
    js = {}
    js['Variables']=[]
    js['Distributions']=[]

    for k,x in enumerate(p):
      js['Variables'].append({})
      js['Variables'][k]['Name'] = x
      js['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

    self.nParameters = len(p)

    k=0
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for beta'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.01
    js['Distributions'][k]['Maximum'] = 5.

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for gamma'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.
    js['Distributions'][k]['Maximum'] = 1.
 
    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for tact'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1.
    js['Distributions'][k]['Maximum'] = 100.
 
    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for dtact'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.0
    js['Distributions'][k]['Maximum'] = 30.0

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for kbeta'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.0
    js['Distributions'][k]['Maximum'] = 1.0



    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for Sigma'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1e-6
    js['Distributions'][k]['Maximum'] = 100.

    return js




  def computational_model( self, s ):
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt,N=N,p=p)

    # get incidents
    y = -np.diff(sol.y[0])
     
    eps = 1e-32
    y[y < eps] = eps
    
    # Transform gradients
    if(self.sampler == 'mTMCMC'):
        sgrad    = []
        diffgrad = []
        for idx in range(len(y)):
            tmp = -(sol.gradMu[idx+1]-sol.gradMu[idx])
            diffgrad.append(list(tmp))
            
            tmp = tmp*p[-1]
            tmp[-1] = y[idx]*sol.gradSig[idx][-1]
            sgrad.append(list(tmp))

        s['Gradient Mean'] = diffgrad
        s['Gradient Standard Deviation'] = sgrad

    s['Reference Evaluations'] = list(y)
    s['Standard Deviation'] = ( p[-1] * y ).tolist()


  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1],t_eval=t.tolist(), N=N,p=p)
    
    y = -np.diff(sol.y[0])
    y = np.append(0, y)

    eps = 1e-32
    y[y < eps] = eps
    
    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][0]['Name']   = 'Daily Incidence'
    js['Variables'][0]['Values'] = list(y)

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    js['Standard Deviation'] = ( p[-1] * y ).tolist()
    s['Saved Results'] = js
