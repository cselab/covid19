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
    
    seir_int  = libepidemics.country.seir_intexp_reparam
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seir_int.Solver(data)

    params = seir_int.Parameters(R0=p[0], D=p[1], Z=p[2], tact=p[3], k=p[4])
    
    s0, i0  = y0
    y0cpp   = (s0, 0.0, i0, 0.0)
    initial = seir_int.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    gradmu    = []
    gradsig   = []

    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S().val()-entry.E().val()
        gradmu.append(np.array([ -entry.S().d(0)-entry.E().d(0), -entry.S().d(1)-entry.E().d(1), -entry.S().d(2)-entry.E().d(2), -entry.S().d(3)-entry.E().d(3), -entry.S().d(4)-entry.E().d(4), -entry.S().d(5)-entry.E().d(5), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = infected
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol

  def computational_model( self, s ):
    p  = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt,N=N,p=p)

    # get incidents
    y   = np.diff(sol.y)
    eps = 1e-32
    y[y < eps] = eps

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
        s['Dispersion'] = ( p[-1] * y ).tolist()
