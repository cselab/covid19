import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.seiir.nbin'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )
 

  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('beta', 0.1, 100), 
            ('mu', 0., 0.1), 
            ('alpha', 0., 0.1),
            ('Z', 0, 50), 
            ('D', 1000, 8000), 
            ('[r]', 0.001, 10)
            )
    
    return js


  def computational_model( self, s ):
 
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt, N=N, p=p)

    # get incidents
    y = -np.diff(sol.y[0])
     
    eps = 1e-32
    y[y < eps] = eps
  
    if(self.sampler == 'mTMCMC'):
        print("[Epidemics] mTMCMC not yet available for nbin")
        sys.exit(0)

    s['Reference Evaluations'] = list(y)
    s['Dispersion'] = ( p[-1] * y ).tolist()


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
    js['Variables'] = [{}]

    js['Variables'][0]['Name'] = 'Daily Incidence'
    js['Variables'][0]['Values'] = list(y)

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    js['Dispersion'] = ( p[-1] * y ).tolist()

    s['Saved Results'] = js
