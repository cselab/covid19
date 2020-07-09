import os
import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):


  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    seird_int_reparam  = libepidemics.country.seird_int_reparam
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seird_int_reparam.Solver(data)

    params = seird_int_reparam.Parameters(R0=p[0], D=p[1], Z=p[2],eps=p[3], tact=p[4], dtact=p[5], kbeta=p[6])
    
    s0, i0  = y0
    y0cpp   = (s0, 0.0, i0, 0.0, 0.0)
    initial = seird_int_reparam.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    exposed   = np.zeros(len(cpp_res))
    deaths   = np.zeros(len(cpp_res))

    gradmu    = []
    gradsig   = []

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S().val()-entry.E().val()
        exposed[idx]   = N-entry.S().val()
        recovered[idx] = entry.R().val()
        deaths[idx] = entry.D().val()

        gradmu.append(np.array([ -entry.S().d(0)-entry.E().d(0), -entry.S().d(1)-entry.E().d(1), -entry.S().d(2)-entry.E().d(2), -entry.S().d(3)-entry.E().d(3), -entry.S().d(4)-entry.E().d(4), -entry.S().d(5)-entry.E().d(5), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.i       = infected
    sol.e       = exposed
    sol.r       = recovered
    sol.d       = deaths
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol

  def computational_model( self, s ):
    p  = s['Parameters']
    t  = self.data['Model']['x-eval']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt,N=N,p=p)

    # get incidences
    infected   = np.diff(sol.i) 
    eps = 1e-32
    infected[infected < eps] = eps
    infected   = infected[self.data['Model']['x-infected']-1] 
    
    # get deaths
    deaths   = np.diff(sol.d)
    eps = 1e-32
    deaths[deaths < eps] = eps
    deaths   = deaths[self.data['Model']['x-deaths']-1] 

    # Concat 
    y = np.concatenate([infected,deaths])
    
    # Transform gradients
    if(self.sampler == 'mTMCMC' and self.likelihoodModel != 'Negative Binomial' ):
        sgrad    = []
        diffgrad = []
        for idx in range(len(y)):
            tmp = (sol.gradMu[idx+1]-sol.gradMu[idx])
            diffgrad.append(list(tmp))
            
            tmp = tmp*p[-1]
            tmp[-1] = y[idx]*sol.gradSig[idx][-1]
            sgrad.append(list(tmp))

        s['Gradient Mean'] = diffgrad
        s['Gradient Standard Deviation'] = sgrad

    elif (self.sampler == 'mTMCMC' and self.likelihoodModel == 'Negative Binomial' ) : 
        print("[Epidemics] mTMCMC not yet available for nbin")
        sys.exit(0)

    s['Reference Evaluations'] = list(y)
    
    if self.likelihoodModel == 'Normal':
        s['Standard Deviation'] = ( p[-1] * y ).tolist()
    elif self.likelihoodModel == 'Positive Normal':
        s['Standard Deviation'] = ( p[-1] * y ).tolist()
    elif self.likelihoodModel == 'Negative Binomial':
        s['Dispersion'] = [p[-1]] * len(y)
