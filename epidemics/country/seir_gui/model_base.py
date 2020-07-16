import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics  #cpp backend


class Object(object):
    pass


class ModelBase(EpidemicsCountry):
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

    def solve_ode(self, y0, T, t_eval, N, p):
 
        seir_int  = libepidemics.country.seir_int
        dp        = libepidemics.country.DesignParameters(N=N)
        cppsolver = seir_int.Solver(dp)
        
        beta = p[0]/5.2
        params = seir_int.Parameters(beta=beta, 
                                     gamma=1/5.2, 
                                     a=1./2.9, 
                                     tact=p[1]-5.0, 
                                     dtact=p[2], 
                                     kbeta=p[3])

        s0, i0 = y0
        y0cpp = (s0, 0.0, i0, 0.0)
        initial = seir_int.State(y0cpp)

        cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt=0.01)

        yS = np.zeros(len(cpp_res))
        gradmu = []
        gradsig = []

        for idx, entry in enumerate(cpp_res):
            yS[idx] = entry.S().val()

        # Fix bad values
        yS[np.isnan(yS)] = 0

        # Create Solution Object
        sol = Object()
        sol.y = yS
        sol.gradMu = gradmu
        sol.gradSig = gradsig

        return sol
