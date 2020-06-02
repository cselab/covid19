import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):


  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    params = sir_int_r0.Parameters(r0=p[0], gamma=1.0/5.2, tact=p[1], dtact=p[2], kbeta=p[3])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    
    # w.r.t. r0, gamma, tact, dtact, kbeta, S0, I0, R0
    params_derivatives = ( 
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0]
    )
    y0_derivatives = (
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    )

    sir_int_r0 = libepidemics.country.sir_int_r0
    
    data   = libepidemics.country.ModelData(N=N)
    solver = sir_int_r0.Solver(data)

    cpp_results, cpp_der_results = country_custom_derivatives(
            solver, params, y0cpp, params_derivatives, y0_derivatives, t_eval=t_eval, dt=0.01)
    
    
    yS      = np.zeros(len(cpp_res))
    gradmu  = []
    gradsig = []

    idx = 0
    for cpp, cpp_der in zip(cpp_results, cpp_der_results):
        yS[idx] = cpp[0]
        gradmu.append(np.array([ cpp_der[0, 0], cpp_der[0, 2], cpp_det[0, 3], cpp_der[0, 4], cpp_der[0, 5], 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))

    # Fix bad values
    yS[np.isnan(yS)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = [yS]
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol
