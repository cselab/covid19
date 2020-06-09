import numpy as np

from epidemics.country.country import EpidemicsCountry

import libepidemics  #cpp backend


class Object(object):
    pass


class ModelBase(EpidemicsCountry):
    def __init__(self, **kwargs):

        super().__init__(**kwargs)

    def solve_ode(self, y0, T, t_eval, N, p):

        sir_int_r0 = libepidemics.country.sir_int_r0
        data = libepidemics.country.ModelData(N=N)
        cppsolver = sir_int_r0.Solver(data)

        params = sir_int_r0.Parameters(r0=p[0],
                                       gamma=1. / 5.2,
                                       tact=p[1],
                                       dtact=p[2],
                                       kbeta=p[3])

        s0, i0 = y0
        y0cpp = (s0, i0, 0.0)
        initial = sir_int_r0.State(y0cpp)

        cpp_res = cppsolver.solve_ad(params, initial, t_eval=t_eval, dt=0.01)

        yS = np.zeros(len(cpp_res))
        gradmu = []
        gradsig = []

        for idx, entry in enumerate(cpp_res):
            yS[idx] = entry.S().val()
            gradmu.append(
                np.array([
                    entry.S().d(0),
                    entry.S().d(1),
                    entry.S().d(2),
                    entry.S().d(3), 0.0
                ]))
            gradsig.append(np.array([0.0, 0.0, 0.0, 0.0, 1.0]))

        # Fix bad values
        yS[np.isnan(yS)] = 0

        # Create Solution Object
        sol = Object()
        sol.y = yS
        sol.gradMu = gradmu
        sol.gradSig = gradsig

        return sol
