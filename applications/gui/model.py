# Created by Petr Karnakov on 25.05.2020
# Copyright 2020 ETH Zurich

import os
import numpy as np

from epidemics.epidemics import EpidemicsBase
from epidemics.tools.tools import save_file
import libepidemics


class Solution():
    def __init__(self, y, gradMu, gradSig):
        self.y = y
        self.gradMu = gradMu
        self.gradSig = gradSig


class Model(EpidemicsBase):
    def __init__(self,
                 dataTotalInfected,
                 populationSize,
                 dataDays=[],
                 params_to_infer=[],
                 params_prior=None,
                 params_fixed=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.dataTotalInfected = dataTotalInfected
        self.dataDays = dataDays if dataDays else np.arange(
            len(dataTotalInfected), dtype=float)
        assert len(self.dataTotalInfected) == len(self.dataDays),\
                "size mismatch {:} != {:}".format(
                        len(self.dataTotalInfected), len(self.dataDays))
        self.populationSize = populationSize
        self.__process_data()

        def update_known(orig, new):
            if not new:
                return
            for k in new:
                assert k in orig, f"Uknown parameter {k}={new[k]}"
                orig[k] = new[k]

        self.params_fixed = {
            'R0': 1.3,
            'I0': 1,
            'gamma': 1. / 5.2,
            'tint': 23.,
            'dint': 10,
            'kint': 0.5,
            '[r]': 1.5,
        }
        update_known(self.params_fixed, params_fixed)

        self.params_prior = {
            'R0': (1., 4.),
            'I0': (0.1, 10.),
            'gamma': (0.09, 0.18),
            'tint': (0., 80.),
            'dint': (1., 30),
            'kint': (0.01, 0.99),
            '[r]': (0., 50.),
        }
        update_known(self.params_prior, params_prior)

        self.likelihoodModel = "Negative Binomial"
        self.params_to_infer = params_to_infer + ['[r]']

        for k in self.params_to_infer:
            assert k in self.params_prior, \
                    "Unknown prior for parameter '{:}'".format(k)

    def save_data_path(self):
        return (self.dataFolder, )

    def get_variables_and_distributions(self):
        self.nParameters = len(self.params_to_infer)
        triples = []
        for k in self.params_to_infer:
            p = self.params_prior[k]
            triples.append((k, *p))
        js = self.get_uniform_priors(*triples)
        return js

    def get_params_to_infer(self):
        return self.params_to_infer

    def __process_data(self):
        t = self.dataDays
        y = self.dataTotalInfected
        N = self.populationSize
        I0 = y[0]
        S0 = N - I0
        y0 = S0, I0

        # FIXME: skip points with daily<=0
        self.data['Model']['x-data'] = t[1:]
        self.data['Model']['y-data'] = np.diff(y)

        self.data['Model']['Initial Condition'] = y0

        save_file(self.data, self.saveInfo['inference data'],
                  'Data for Inference', 'pickle')

    def sample(self, nSamples):
        super().sample(nSamples, 0.7)

    def propagate(self, nPropagate, futureDays):
        T = np.ceil(self.data['Model']['x-data'][-1] + futureDays)
        self.data['Propagation']['x-data'] = np.linspace(0., T, int(T + 1))
        super().propagate(nPropagate)

    def substitute_inferred(self, p):
        """
        p: `list(float)`
            Parameter values from Korali
        Returns:
        `dict(str -> float)`
            Inferred and fixed parameters.
        """

        r = dict(self.params_fixed)
        r.update({var: p[i] for i, var in enumerate(self.params_to_infer)})
        return r

    def __solve_ode(self, y0, T, t_eval, N, p):

        sir_int_r0 = libepidemics.country.sir_int_r0
        dp = libepidemics.country.DesignParameters(N=N)
        cppsolver = sir_int_r0.Solver(dp)

        pp = self.substitute_inferred(p)

        params = sir_int_r0.Parameters(r0=pp['R0'],
                                       gamma=pp['gamma'],
                                       tact=pp['tint'],
                                       dtact=pp['dint'],
                                       kbeta=pp['kint'])

        s0, i0 = y0
        y0cpp = (s0, i0, 0.0)
        initial = sir_int_r0.State(y0cpp)

        AUTOGRAD = False
        if AUTOGRAD:
            cpp_res = cppsolver.solve_params_ad(params,
                                                initial,
                                                t_eval=t_eval,
                                                dt=0.1)

            yS = np.zeros(len(cpp_res))
            gradmu = []  # gradients wrt model parameters
            gradsig = []  # gradients wrt distribution parameter

            for i, entry in enumerate(cpp_res):
                yS[i] = entry.S().val()
                gradmu.append(
                    np.array([
                        entry.S().d(0),
                        entry.S().d(1),
                        entry.S().d(2),
                        entry.S().d(3), 0.0
                    ]))
                gradsig.append(np.array([0.0, 0.0, 0.0, 0.0, 1.0]))
        else:
            cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt=0.1)

            yS = np.zeros(len(cpp_res))
            gradmu = []  # gradients wrt model parameters
            gradsig = []  # gradients wrt distribution parameter

            for i, entry in enumerate(cpp_res):
                yS[i] = entry.S()

        return Solution(y=[yS], gradMu=gradmu, gradSig=gradsig)

    def computational_model(self, s):
        p = s['Parameters']
        t = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.populationSize

        tt = [t[0] - 1] + t.tolist()
        sol = self.__solve_ode(y0=y0, T=t[-1], t_eval=tt, N=N, p=p)

        # daily
        y = -np.diff(sol.y[0])

        eps = 1e-32
        y[y < eps] = eps

        if (self.sampler == 'mTMCMC'):
            raise RuntimeError("mTMCMC not yet available for nbin")

        s['Reference Evaluations'] = list(y)
        s['Dispersion'] = [p[-1]] * len(y)

    def computational_model_propagate(self, s):
        p = s['Parameters']
        t = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.populationSize

        tt = [t[0] - 1] + t.tolist()
        sol = self.__solve_ode(y0=y0, T=t[-1], t_eval=tt, N=N, p=p)

        y = -np.diff(sol.y[0])

        eps = 1e-32
        y[y < eps] = eps

        js = {}
        js['Variables'] = []

        js['Variables'].append({})
        js['Variables'][0]['Name'] = 'Daily Incidence'
        js['Variables'][0]['Values'] = list(y)

        js['Number of Variables'] = len(js['Variables'])
        js['Length of Variables'] = len(t)

        js['Dispersion'] = [p[-1]] * len(y)

        s['Saved Results'] = js
