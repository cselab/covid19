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
    
    seird_int_reparam  = libepidemics.country.seird_intsmooth_reparam
    dp        = libepidemics.country.DesignParameters(N=N)
    cppsolver = seird_int_reparam.Solver(dp)

    params = seird_int_reparam.Parameters(R0=p[0], D=1.0/self.constants['gamma'], Z=p[1], eps=p[2], tact=p[3], dtact=p[4], kbeta=p[5])
    
    s0, i0  = y0
    y0cpp   = (s0, 0.0, i0, 0.0, 0.0) # S E I R D
    initial = seird_int_reparam.State(y0cpp)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    exposed   = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S()-entry.E()
        exposed[idx]   = N-entry.S()
        recovered[idx] = entry.R()
        deaths[idx]    = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = infected
    sol.e       = exposed
    sol.r       = recovered
    sol.d       = deaths
 
    return sol
