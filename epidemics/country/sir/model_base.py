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
    
    sir       = libepidemics.country.sir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = sir.Solver(data)

    params = sir.Parameters(beta=p[0], gamma=p[1])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    initial = sir.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.1)
    
    yS      = np.zeros(len(cpp_res))
    gradmu  = []
    gradsig = []

    for idx,entry in enumerate(cpp_res):
        yS[idx] = entry.S().val()
        gradmu.append(np.array([ entry.S().d(0), entry.S().d(1), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 1.0 ]))

    # Fix bad values
    yS[np.isnan(yS)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = yS
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol
