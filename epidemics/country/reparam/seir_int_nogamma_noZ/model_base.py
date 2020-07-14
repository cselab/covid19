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
    
    seir_int  = libepidemics.country.seir_int_reparam
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seir_int.Solver(data)

    params = seir_int.Parameters(R0=p[0], D=1.0/self.constants['gamma'], Z=self.constants['Z'], tact=p[1], dtact=p[2], kbeta=p[3])
    
    s0, i0  = y0
    y0cpp   = (s0, 0.0, i0, 0.0)
    initial = seir_int.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected  = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    exposed   = np.zeros(len(cpp_res))
    gradmu    = []
    gradsig   = []

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S().val()-entry.E().val()
        exposed[idx]   = N-entry.S().val()
        recovered[idx] = entry.R().val()
        gradmu.append(np.array([ -entry.S().d(0)-entry.E().d(0), -entry.S().d(1)-entry.E().d(1), -entry.S().d(2)-entry.E().d(2), -entry.S().d(3)-entry.E().d(3), -entry.S().d(4)-entry.E().d(4), -entry.S().d(5)-entry.E().d(5), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))

    # Fix bad values
    infected[np.isnan(infected)] = 0
    
    # Create Solution Object
    sol = Object()
    sol.y       = infected
    sol.e       = exposed
    sol.r       = recovered
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol
