# Created by Petr Karnakov on 25.05.2020
# Copyright 2020 ETH Zurich

import os
import numpy as np

from epidemics.epidemics import EpidemicsBase
from epidemics.tools.tools import save_file
import libepidemics

class Object():
    pass


class Model(EpidemicsBase):
    def __init__(self,
                 data,
                 populationSize,
                 futureDays=0,
                 nSamples=2000,
                 nPropagation=100,
                 **kwargs):
        super().__init__(**kwargs)

        self.futureDays = futureDays
        self.totalInfected = data
        self.populationSize = populationSize
        self.__process_data()
        self.likelihoodModel = "Negative Binomial"

    def save_data_path(self):
        return (self.dataFolder, )

    def get_variables_and_distributions(self):
        self.nParameters = 4
        js = self.get_uniform_priors(
            ('R0', 1, 4),
            ('tact', 0, 80),
            ('kbeta', 0.1, 10),
            ('[r]', 0, 50),
        )
        return js

    def __process_data(self):
        t = np.arange(60)
        y = t + 1
        N = 1e6
        I0 = y[0]
        S0 = N - I0
        y0 = S0, I0

        self.data['Model']['x-data'] = t[1:]
        self.data['Model']['y-data'] = np.diff(y[0:])

        self.data['Model']['Initial Condition'] = y0
        self.data['Model']['Population Size'] = self.populationSize

        T = np.ceil(t[-1] + self.futureDays)
        self.data['Propagation']['x-data'] = np.linspace(0, T, int(T + 1))

        save_file(self.data, self.saveInfo['inference data'],
                  'Data for Inference', 'pickle')

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

        cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt=0.01)

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
        sol.y = [yS]
        sol.gradMu = gradmu
        sol.gradSig = gradsig

        return sol

    def computational_model(self, s):
        p = s['Parameters']
        t = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']

        tt = [t[0] - 1] + t.tolist()
        sol = self.solve_ode(y0=y0, T=t[-1], t_eval=tt, N=N, p=p)
        y = -np.diff(sol.y[0])

        # get incidents
        y = -np.diff(sol.y[0])

        eps = 1e-32
        y[y < eps] = eps

        if (self.sampler == 'mTMCMC'):
            raise RuntimeError("mTMCMC not yet available for nbin")

        s['Reference Evaluations'] = list(y)
        s['Dispersion'] = (p[-1] * y).tolist()

    def computational_model_propagate(self, s):
        p = s['Parameters']
        t = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']

        tt = [t[0] - 1] + t.tolist()
        sol = self.solve_ode(y0=y0, T=t[-1], t_eval=t.tolist(), N=N, p=p)

        y = -np.diff(sol.y[0])
        y = np.append(0, y)

        eps = 1e-32
        y[y < eps] = eps

        js = {}
        js['Variables'] = []

        js['Variables'].append({})
        js['Variables'][0]['Name'] = 'Daily Incidence'
        js['Variables'][0]['Values'] = list(y)

        js['Number of Variables'] = len(js['Variables'])
        js['Length of Variables'] = len(t)

        js['Dispersion'] = (len(y)) * [p[-1]]

        s['Saved Results'] = js
