import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):


  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    sird_int   = libepidemics.country.sird_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = sird_int.Solver(dp)

    params = sird_int.Parameters(R0=p[0], D=p[1], eps=p[2], tact=self.intday+p[3], dtact=p[4], kbeta=p[5])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0, 0.0) # S I R D
    initial = sird_int.State(y0cpp)
    
    cpp_res = None
    if(self.sampler == 'mTMCMC' or self.sampler=='HMC'):
        cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.1)
    else:
        cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.1)
    
    infected   = np.zeros(len(cpp_res))
    recovered  = np.zeros(len(cpp_res))
    deaths     = np.zeros(len(cpp_res))
    gradinfected = []
    graddeaths = []
    gradsig    = []

    for idx,entry in enumerate(cpp_res):
        if(self.sampler == 'mTMCMC' or self.sampler=='HMC'):
            infected[idx]  = N-entry.S().val()
            recovered[idx] = entry.R().val()
            deaths[idx]    = entry.D().val()
            gradinfected.append(np.array([ -entry.S().d(0), -entry.S().d(1), 0.0 ])) 
            graddeaths.append(np.array([ entry.D().d(0), entry.D().d(1), 0.0 ])) 
            gradsig.append(np.array([ 0.0, 0.0, 1.0 ]))

        else:
            infected[idx]  = N-entry.S()
            recovered[idx] = entry.R()
            deaths[idx]    = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0
    
    # Create Solution Object
    sol = Object()
    sol.y = infected
    sol.r = recovered
    sol.d = deaths
    sol.gradMu  = gradinfected
    sol.gradSig = gradsig
 
    return sol
