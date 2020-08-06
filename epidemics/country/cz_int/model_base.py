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
    
    cz_int    = libepidemics.country.cz_int
    dp        = libepidemics.country.DesignParameters(N=N)
    cppsolver = cz_int.Solver(dp)

    D = p[1]
    Z = p[2]

    params = cz_int.Parameters(R0=p[0], 
            gamma=1/D, sigma=1/Z,
            eps1=self.bz_constants['eps1'], eps2=self.bz_constants['eps2'],
            eps3=p[3], eps4=self.bz_constants['eps4'],
            omega1=self.bz_constants['omega1'], omega2=self.bz_constants['omega2'],
            omega3=self.bz_constants['omega3'], omega4=self.bz_constants['omega4'],
            omega5=self.bz_constants['omega5'], 
            tact=self.intday+p[4], dtact=p[5], kbeta=p[6]) # 7 params
    
    s0, i0  = y0
    y0cpp   = (s0, p[0]*i0, i0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # S E I P H1 H2 U R D C
    
    initial = cz_int.State(y0cpp)
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    exposed   = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = entry.C()
        exposed[idx]   = N-entry.S()
        recovered[idx] = entry.R()
        deaths[idx]    = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = infected
    sol.e       = exposed
    sol.r       = recovered
    sol.d       = deaths
 
    return sol
