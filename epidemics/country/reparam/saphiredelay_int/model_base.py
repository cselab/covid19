import math
import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsCountry ):


  def __init__( self, **kwargs ):

    super().__init__( **kwargs )

  def solve_ode( self, y0, T, t_eval, N, p ):
    
    saphire_int = libepidemics.country.saphire_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = saphire_int.Solver(dp)

    params = saphire_int.Parameters(R0=p[0], D=p[1], Z=p[2], Y=p[3], mu=p[4], alpha=p[5], eps=p[6], tact=self.intday+p[7], dtact=p[8], kbeta=p[9])
  
    beta = p[0]/p[1]

    s0, ir0  = y0
    iu0 = (1-p[5])/p[5]*ir0
    i0  = iu0 + ir0
    
    p0 = beta*p[3]*i0
    e0 = beta*p[2]*p0
    s0 = s0 - e0 - p0 - i0

    y0cpp   = (s0, e0, p0, ir0, iu0, 0.0, 0.0, 0.0, 0.0) # S E P Ir Iu R D Cir Ciu
    
    initial = saphire_int.State(y0cpp)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected        = np.zeros(len(cpp_res))
    exposed         = np.zeros(len(cpp_res))
    infectedu       = np.zeros(len(cpp_res))
    recovered       = np.zeros(len(cpp_res))
    preasymptomatic = np.zeros(len(cpp_res))
    deaths          = np.zeros(len(cpp_res))
    
    cir = np.zeros(len(cpp_res))
    ciu = np.zeros(len(cpp_res))

    dt = p[10]
    w1 = math.ceil(dt)-dt
    w2 = 1.-w1


    for idx,entry in enumerate(cpp_res):
        exposed[idx]         = N-entry.S()
        preasymptomatic[idx] = N-entry.S()-entry.E()
        infected[idx]        = (N-entry.S()-entry.E()-entry.P())*p[5]
        infectedu[idx]       = (N-entry.S()-entry.E()-entry.P())*(1-p[5])
        recovered[idx]       = entry.R()
        
        cir[idx] = entry.Cir()
        ciu[idx] = entry.Ciu()
        
        if math.floor(idx+dt) < len(deaths):
            deaths[math.floor(idx+dt)] += w1 * entry.D()
        if math.ceil(idx+dt) < len(deaths):
            deaths[math.ceil(idx+dt)]  += w2 * entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0
    
    # Create Solution Object
    sol = Object()
    sol.y = infected
    sol.iu = infectedu
    sol.p = preasymptomatic
    sol.e = exposed
    sol.r = recovered
    sol.d = deaths
    sol.cir = cir
    sol.ciu = ciu
 
    return sol
