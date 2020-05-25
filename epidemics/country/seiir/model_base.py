import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):

  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):

    seiir     = libepidemics.country.seiir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seiir.Solver(data)

    params = seiir.Parameters(beta=p[0], mu=p[1], alpha=p[2], Z=p[3], D=p[4])
 
    s0, ir0 = y0
    y0cpp   = (s0, 0.0, ir0, 0.0, 0.0) # S E Ir Iu  R
    
    initial = seiir.State(y0cpp)
 
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
  
    yS      = np.zeros(len(cpp_res))
    gradmu  = []
    gradsig = []

    for idx,entry in enumerate(cpp_res):
        yS[idx] = entry.S().val()
        gradmu.append(np.array([ entry.S().d(0), entry.S().d(1), entry.S().d(2), entry.S().d(3), entry.S().d(4), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))
 
    # Fix bad values
    yS[np.isnan(yS)] = 0
 
    # Create Solution Object
    sol = Object()
    sol.y       = [yS]
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol
