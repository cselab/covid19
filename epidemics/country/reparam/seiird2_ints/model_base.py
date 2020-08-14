import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):

  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    seiird2_ints = libepidemics.country.seiird2_ints_reparam
    dp           = libepidemics.country.DesignParameters(N=N)
    cppsolver    = seiird2_ints.Solver(dp)

    params = seiird2_ints.Parameters(R0=p[0], D=p[1], Z=p[2], mu=p[3], alpha=p[4], eps=p[5], tact=self.intday+p[6], kbeta=p[7])

    s0, ir0  = y0
    iu0  = (1-p[4])/p[4]*ir0
    e0   = (np.exp(p[0])-1)/np.exp(1/p[2]) * ir0 + (np.exp(p[3]*p[0])-1)/np.exp(1/p[2]) * iu0

    y0cpp   = (s0, e0, ir0, iu0, 0.0, 0.0) # S E Ir Iu  R D
   
    initial = seiird2_ints.State(y0cpp)
 
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
  
    exposed   = np.zeros(len(cpp_res))
    infected  = np.zeros(len(cpp_res))
    infectedu = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        exposed[idx]   = N-entry.S()
        infected[idx]  = N-entry.S()-entry.E()-entry.Iu()
        infectedu[idx] = N-entry.S()-entry.E()-entry.Ir()
        recovered[idx] = entry.R()
        deaths[idx]    = entry.D()
        
 
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
