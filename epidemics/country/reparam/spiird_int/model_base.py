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
    
    spiird_int = libepidemics.country.spiird_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = spiird_int.Solver(dp)

    params = spiird_int.Parameters(R0=p[0], D=p[1], Y=p[2], alpha=p[3], eps=p[4], tact=self.intday+p[5], dtact=p[6], kbeta=p[7])
    
    s0, ir0 = y0
    iu0     = (1-p[3])/p[3] * ir0
    e0      = (np.exp(p[0])-1)/np.exp(1/p[2]) * (iu0 + Ir0)
 
    y0cpp   = (s0, e0, ir0, iu0, 0.0, 0.0) # S P Ir Iu R D
    
    initial = spiird_int.State(y0cpp)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected        = np.zeros(len(cpp_res))
    infectedu       = np.zeros(len(cpp_res))
    recovered       = np.zeros(len(cpp_res))
    preasymptomatic = np.zeros(len(cpp_res))
    deaths          = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        preasymptomatic[idx] = N-entry.S()
        infected[idx]        = N-entry.S()-entry.P()-entry.Iu()
        infectedu[idx]       = N-entry.S()-entry.P()-entry.Ir()
        recovered[idx]       = entry.R()
        deaths[idx]          = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0
    
    # Create Solution Object
    sol = Object()
    sol.y = infected
    sol.p = preasymptomatic
    sol.r = recovered
    sol.d = deaths
 
    return sol
