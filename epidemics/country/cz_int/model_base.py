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

    params = cz_int.Parameters(R0=p[0], 
            gamma=self.bz_constants['gamma'], sigma=self.bz_constants['sigma'],
            eps1=self.bz_constants['eps1'], eps2=self.bz_constants['eps2'],
            eps3=self.bz_constants['eps3'], eps4=self.bz_constants['eps4'],
            omega1=self.bz_constants['omega1'], omega2=self.bz_constants['omega2'],
            omega3=self.bz_constants['omega3'], omega4=self.bz_constants['omega4'],
            omega5=self.bz_constants['omega5'], 
            tact=p[1], dtact=p[2], kbeta=p[3])
    
    s0, i0  = y0
    y0cpp   = (s0, 0.0, i0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) # S E I P H1 H2 U R D C
    
    initial = cz_int.State(y0cpp)
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    exposed   = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = entry.C().val()
        exposed[idx]   = N-entry.S().val()
        recovered[idx] = entry.R().val()
        deaths[idx]    = entry.D().val()

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
