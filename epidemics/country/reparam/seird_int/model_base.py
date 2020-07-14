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
    y0cpp   = (s0, 0.0, i0, 0.0, 0.0) # S E I R D
    initial = seird_int_reparam.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    exposed   = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S().val()-entry.E().val()
        exposed[idx]   = N-entry.S().val()
        recovered[idx] = entry.R().val()
        deaths[idx]    = entry.D().val()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = infected
    sol.e       = exposed
    sol.r       = recovered
    sol.d       = deaths
 
    return sol