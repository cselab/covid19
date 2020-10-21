import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):

  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    seiird2_int = libepidemics.country.seiird2_int_reparam
    dp          = libepidemics.country.DesignParameters(N=N)
    cppsolver   = seiird2_int.Solver(dp)

    params = seiird2_int.Parameters(R0=p[0], D=p[1], Z=p[2], mu=p[3], alpha=p[4], eps=p[5], tact=self.intday+p[6], dtact=p[7], kbeta=p[8])
 
    beta = p[0]/p[1]
    
    s0, ir0  = y0
    iu0 = (1-p[4])/p[4]*ir0
    i0  = ir0 + iu0

    e0 = beta*p[2]*i0
    s0 = s0 - e0 - iu0
 
    y0cpp  = (s0, e0, ir0, iu0, 0.0, 0.0) # S E Ir Iu  R D
    
    initial = seiird2_int.State(y0cpp)
 
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
  
    exposed   = np.zeros(len(cpp_res))
    infected  = np.zeros(len(cpp_res))
    infectedu = np.zeros(len(cpp_res))
    recovered = np.zeros(len(cpp_res))
    deaths    = np.zeros(len(cpp_res))

    cir = np.zeros(len(cpp_res))
    ciu = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        exposed[idx]   = N-entry.S()
        infected[idx]  = N-entry.S()-entry.E()-entry.Iu()
        infectedu[idx] = N-entry.S()-entry.E()-entry.Ir()
        recovered[idx] = entry.R()
        deaths[idx]    = entry.D()
        cir[idx] = entry.Cir()
        ciu[idx] = entry.iu()
        
 
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0
 
    # Create Solution Object
    sol = Object()
    sol.y  = infected
    sol.iu = infectedu
    sol.e  = exposed
    sol.r  = recovered
    sol.d  = deaths
    sol.cir = cir
    sol.ciu = ciu
 
    return sol
