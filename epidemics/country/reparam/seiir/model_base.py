import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):

  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):

    seiir     = libepidemics.country.seiir_reparam
    dp        = libepidemics.country.DesignParameters(N=N)
    cppsolver = seiir.Solver(dp)

    params = seiir.Parameters(R0=p[0], D=p[1], Z=p[2],mu=p[3], alpha=p[4])
 
    s0, ir0 = y0
    y0cpp   = (s0, 0.0, ir0, 0.0, 0.0) # S E Ir Iu  R
    
    initial = seiir.State(y0cpp)
 
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.1)
  
    infected = np.zeros(len(cpp_res))
    gradmu   = []
    gradsig  = []

    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S().val()-entry.E().val()-entry.Iu().val()
        gradmu.append(np.array([ 
            -entry.S().d(0)-entry.E().d(0)-entry.Iu().d(0),
            -entry.S().d(1)-entry.E().d(1)-entry.Iu().d(1),
            -entry.S().d(2)-entry.E().d(2)-entry.Iu().d(2),
            -entry.S().d(3)-entry.E().d(3)-entry.Iu().d(3),
            -entry.S().d(4)-entry.E().d(4)-entry.Iu().d(4),
            0.0]))

        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))
 
    infected[np.isnan(infected)] = 0
 
    # Create Solution Object
    sol = Object()
    sol.y       = infected
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol
