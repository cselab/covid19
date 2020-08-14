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
    
    seirud_ints = libepidemics.country.seirud_ints_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = seirud_ints.Solver(dp)

    params = seirud_ints.Parameters(R0=p[0], D=p[1], Z=p[2], Y=p[3], alpha=p[4], eps=self.constants['eps'], tact=self.intday+p[5], kbeta=p[6])
 
    s0, ir0 = y0
    iu0     = (1-p[4])/p[4] * ir0
    p0      = (np.exp(p[0])-1)/np.exp(1/p[3]) * iu0
    e0      = (np.exp(p[0])-1)/np.exp(1/p[2]) * p0

    y0cpp   = (s0, e0, p0, ir0, iu0, 0.0, 0.0) # S E P Ir Iu R D
    
    initial = seirud_int.State(y0cpp)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected        = np.zeros(len(cpp_res))
    exposed         = np.zeros(len(cpp_res))
    infectedu       = np.zeros(len(cpp_res))
    recovered       = np.zeros(len(cpp_res))
    preasymptomatic = np.zeros(len(cpp_res))
    deaths          = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        exposed[idx]         = N-entry.S()
        preasymptomatic[idx] = N-entry.S()-entry.E()
        infected[idx]        = N-entry.S()-entry.E()-entry.P()-entry.Iu()
        infectedu[idx]       = N-entry.S()-entry.E()-entry.P()-entry.Ir()
        recovered[idx]       = entry.R()
        deaths[idx]          = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0
    
    # Create Solution Object
    sol = Object()
    sol.y = infected
    sol.p = preasymptomatic
    sol.e = exposed
    sol.r = recovered
    sol.d = deaths
 
    return sol
