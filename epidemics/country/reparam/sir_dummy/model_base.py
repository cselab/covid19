import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):


  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    sir_int   = libepidemics.country.sir_int_reparam
    dp        = libepidemics.country.DesignParameters(N=N)
    cppsolver = sir_int.Solver(dp)

    params = sir_int.Parameters(R0=p[0], D=p[1], tact=p[2], dtact=p[3], kbeta=0.3)
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    initial = sir_int.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.1)
    
    infected   = np.zeros(len(cpp_res))
    recovered  = np.zeros(len(cpp_res))
    gradmu  = []
    gradsig = []
    graddisp = []

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S().val()
        recovered[idx] = entry.R().val()
        gradmu.append(np.array([ -entry.S().d(0), -entry.S().d(1), -entry.S().d(2), -entry.S().d(3),  0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, infected[idx] ]))
        graddisp.append(np.array([ 0.0, 0.0, 0.0, 0.0, 1.0 ]))

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y        = infected
    sol.r        = recovered
    sol.gradMu   = gradmu
    sol.gradSig  = gradsig
    sol.gradDisp = graddisp
 
    return sol
