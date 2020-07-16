import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):

  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):

    seiird_int = libepidemics.country.seiird_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = seiird_int.Solver(dp)

    params = seiird_int.Parameters(R0=p[0], D=1.0/self.constants['gamma'], Z=p[1], mu=p[2], alpha=p[3], eps=p[4], tact=p[5], dtact=p[6], kbeta=p[7])

    s0, ir0 = y0
    y0cpp   = (s0, 0.0, ir0, 0.0, 0.0, 0.0) # S E Ir Iu R D
    
    initial = seiird_int.State(y0cpp)
 
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
  
    exposed   = np.zeros(len(cpp_res))
    infected  = np.zeros(len(cpp_res))
    infectedu = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        exposed[idx]   = N-entry.S().val()
        infected[idx]  = N-entry.S().val()-entry.E().val()-entry.Iu().val()
        infectedu[idx] = N-entry.S().val()-entry.E().val()-entry.Ir().val()
        recovered[idx] = entry.R().val()
        deaths[idx]    = entry.D().val()
        
 
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0
 
    # Create Solution Object
    sol = Object()
    sol.y  = infected
    sol.iu = infectedu
    sol.e  = exposed
    sol.r  = recovered
    sol.d  = deaths
 
    return sol
