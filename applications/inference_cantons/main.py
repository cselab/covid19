#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import copy
from epidemics.tools.tools import import_from
import argparse

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.tools.tools import save_file
from epidemics.epidemics import EpidemicsBase

from ode import Ode, Sir, Seir, SeirCpp
from data import Data, get_all_canton_keys, get_data_switzerland,\
        get_data_switzerland_cantons, get_data_synthetic


def repeat(v, numRegions):
    return list(np.array([v] * numRegions).T.flatten())


class Model(EpidemicsBase):
    def __init__(self,
                 data: Data,
                 ode: Ode,
                 params_to_infer: list,
                 nPropagation=100,
                 percentages=[0.5, 0.95, 0.99],
                 logPlot=False,
                 **kwargs):
        self.modelName = 'cantons'
        self.modelDescription = 'Fit SEIR* with cantons on Daily Infected Data'
        self.likelihoodModel = 'Negative Binomial'
        self.nPropagation = nPropagation
        self.percentages = percentages
        self.logPlot = logPlot
        self.ode = ode
        self.params_to_infer = params_to_infer

        self.propagationData = {}

        super().__init__(**kwargs)

        self.process_data(data)

    def save_data_path(self):
        return (self.dataFolder, self.modelName)

    def process_data(self, data: Data):
        Itotal = np.array(data.total_infected)  # shape (n_regions, nt)
        t = np.array(data.time)
        N = np.array(data.population)
        I0 = Itotal[:, 0] * 0  # XXX
        #I0 = Itotal[:, 0] * 0 + 1 # XXX
        S0 = N - I0
        y0 = S0, I0

        self.n_regions = Itotal.shape[0]
        self.region_names = data.name

        Idaily = np.diff(Itotal, axis=1)  # shape (n_regions, nt-1)
        Idaily = Idaily.T  # shape (nt-1, n_regions)

        # XXX adhoc, negative binomial requires positive values
        Idaily = np.maximum(1e-10, Idaily)
        #Idaily = np.maximum(0.1, Idaily)

        self.data['Model']['x-data'] = repeat(t[1:], self.n_regions)
        self.data['Model']['y-data'] = Idaily.flatten()

        self.data['Model']['Initial Condition'] = y0
        self.data['Model']['Population Size'] = N

        if data.fit_importance is not None:
            F = np.array(data.fit_importance)
        else:
            F = np.ones([self.n_regions])
        self.data['Model']['Fit Importance'] = F

        if data.commute_matrix is not None:
            C = np.array(data.commute_matrix)
        else:
            C = np.zeros([self.n_regions] * 2)
        self.data['Model']['Commute Matrix'] = C

        if data.commute_airports is not None:
            Qa = np.array(data.commute_airports)
        else:
            Qa = np.zeros([self.n_regions])
        self.data['Model']['Infected Commuters Airports'] = Qa

        if data.commute_borders is not None:
            Qb = np.array(data.commute_borders)
        else:
            Qb = np.zeros([self.n_regions])
        self.data['Model']['Infected Commuters Borders'] = Qb

        T = np.ceil(t[-1])
        self.data['Propagation']['x-data'] = repeat(
            np.linspace(1, T, int(T + 1)), self.n_regions)

        save_file(self.data, self.saveInfo['inference data'],
                  'Data for Inference', 'pickle')

    def get_variables_and_distributions(self):
        p = self.params_to_infer + ['[r]']
        js = {}
        js['Variables'] = []
        js['Distributions'] = []
        for k, name in enumerate(p):
            js['Variables'].append({})
            js['Variables'][k]['Name'] = name
            js['Variables'][k]['Prior Distribution'] = 'Prior for ' + name

        self.nParameters = len(p)

        k = 0

        for name in self.params_to_infer:
            js['Distributions'].append({})
            js['Distributions'][k]['Name'] = 'Prior for ' + name
            js['Distributions'][k]['Type'] = 'Univariate/Uniform'
            minmax = self.ode.params_prior[name]
            js['Distributions'][k]['Minimum'] = minmax[0]
            js['Distributions'][k]['Maximum'] = minmax[1]
            k += 1

        js['Distributions'].append({})
        js['Distributions'][k]['Name'] = 'Prior for [r]'
        js['Distributions'][k]['Type'] = 'Univariate/Uniform'
        js['Distributions'][k]['Minimum'] = 0.01
        js['Distributions'][k]['Maximum'] = 10.
        k += 1
        return js

    def get_params(self, korali_p, N, C, Qa, Qb):
        """
        N: `array_like`, (n_regions)
        Population of regions.
        C: `array_like`, (n_regions,n_regions)
        Commute matrix.
        """
        params = dict(self.ode.params_fixed.items())
        for i, name in enumerate(self.params_to_infer):
            params[name] = korali_p[i]
        params['N'] = np.array(N)
        params['C'] = np.array(C)
        params['Qa'] = np.array(Qa)
        params['Qb'] = np.array(Qb)
        return params

    def computational_model(self, s):
        p = s['Parameters']
        t = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']
        C = self.data['Model']['Commute Matrix']
        Qa = self.data['Model']['Infected Commuters Airports']
        Qb = self.data['Model']['Infected Commuters Borders']
        F = self.data['Model']['Fit Importance']

        tt1 = [min(t) - 1] + list(t[0::self.n_regions])
        assert min(tt1) >= 0

        params = self.get_params(p, N, C, Qa, Qb)
        S, _ = self.ode.solve_SI(params, [0, max(t)], y0, tt1)
        # S: shape (n_regions, nt)
        Idaily = -np.diff(S, axis=1)  # shape (n_regions, nt-1)
        Idaily = Idaily.T  # shape (nt-1, n_regions)
        Idaily = list(Idaily.flatten())

        s['Reference Evaluations'] = Idaily
        dispersion = np.ones(len(Idaily)) * p[-1]
        for i in range(self.n_regions):
            dispersion[i::self.n_regions] *= F[i]
        s['Dispersion'] = list(dispersion)

    def computational_model_propagate(self, s):
        p = s['Parameters']
        t = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']
        C = self.data['Model']['Commute Matrix']
        Qa = self.data['Model']['Infected Commuters Airports']
        Qb = self.data['Model']['Infected Commuters Borders']
        F = self.data['Model']['Fit Importance']

        t1 = list(t[0::self.n_regions])
        assert min(t1) >= 0

        params = self.get_params(p, N, C, Qa, Qb)
        S, _ = self.ode.solve_SI(params, [0, max(t)], y0, t1)
        # S: shape (n_regions, nt)
        # shape (n_regions, nt-1)
        Idaily = -np.diff(S, axis=1)
        # shape (n_regions, nt)
        Idaily = np.hstack((np.zeros((self.n_regions, 1)), (Idaily)))

        js = {}
        js['Variables'] = []

        # variables that will be visible in self.propagatedVariables
        # (does not affect sampling)
        for i in range(self.n_regions):
            js['Variables'].append({})
            js['Variables'][i]['Name'] = 'Daily Incidence {:}'.format(i)
            js['Variables'][i]['Values'] = list(Idaily[i])

        js['Number of Variables'] = len(js['Variables'])
        js['Length of Variables'] = len(t1)

        dispersion = np.ones(len(t1)) * p[-1]
        for i in range(self.n_regions):
            dispersion[i::self.n_regions] *= F[i]
        js['Dispersion'] = list(dispersion)

        s['Saved Results'] = js

    def evaluate(self, korali_p):
        """
        Evaluates one sample of the model with given set of parametrs.
        korali_p: `list`
            Parameters that would be passed
            to`computational_model_propagate` in `s['Parameters']`.
        Output:
            Fills `self.propagatedVariables` with one sample.
        """
        s = dict()
        s['Parameters'] = korali_p
        self.computational_model_propagate(s)
        js = s['Saved Results']
        self.propagatedVariables = {}
        Ns = 1
        Nv = js['Number of Variables']
        Nt = js['Length of Variables']
        for i in range(len(js['Variables'])):
            varName = js['Variables'][i]['Name']
            self.propagatedVariables[varName] = np.zeros((Ns, Nt))
            for k in range(Ns):
                self.propagatedVariables[varName][k] = np.asarray(
                    js['Variables'][i]['Values'])
        varName = 'Dispersion'
        self.propagatedVariables[varName] = np.zeros((Ns, Nt))
        for k in range(Ns):
            self.propagatedVariables[varName][k] = np.asarray(js['Dispersion'])


def main():
    nSamples = 1000
    x = argparse.Namespace()
    x.dataFolder = "data/"
    x.nPropagation = 20
    x.percentages = [0.5]
    x.nThreads = 8

    #data = get_data_switzerland()
    #data = get_data_synthetic()

    keys = get_all_canton_keys()
    #keys = ['ZH', 'AG', 'TI']

    data = get_data_switzerland_cantons(keys)

    data.fit_importance = [1] * len(keys)
    key_sel = []
    #key_sel = ['ZH', 'TI']
    for k in key_sel:
        data.fit_importance[keys.index(k)] *= 100

    #ode = Sir()
    #params_to_infer = ['R0', 'gamma']

    #ode = Seir()
    ode = SeirCpp()
    #params_to_infer = []
    params_to_infer = ['R0', 'Z', 'D', 'nu', 'theta_b', 'tact']
    ode.params_fixed['theta_a'] = 0.

    #ode.params_fixed['tact'] = 1e10
    #params_to_infer = ['R0', 'Z', 'D', 'nu', 'theta_b']

    #ode.params_fixed['nu'] = 0
    #params_to_infer = ['R0', 'Z', 'D', 'theta_b', 'tact']

    a = Model(data, ode, params_to_infer, **vars(x))

    if 1:
        a.sample(nSamples)
        a.propagate()
    else:
        a.evaluate([1.34, 1., 2., 2.5])

    a.save()


if __name__ == "__main__":
    main()
