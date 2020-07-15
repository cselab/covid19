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
    data       = libepidemics.country.ModelData(N=N)
    cppsolver  = sird_int.Solver(data)

    params = sird_int.Parameters(R0=p[0], D=1.0/self.constants['gamma'], eps=p[1], tact=p[2], dtact=p[3], kbeta=p[4])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0, 0.0) # S I R D
    initial = sird_int.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected   = np.zeros(len(cpp_res))
    recovered  = np.zeros(len(cpp_res))
    deaths     = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S().val()
        recovered[idx] = entry.R().val()
        deaths[idx]    = entry.D().val()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y = infected
    sol.r = recovered
    sol.d = deaths
 
    return sol
