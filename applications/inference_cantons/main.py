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
import epidemics.data.swiss_cantons as swiss_cantons


def repeat(v, numRegions):
    return list(np.array([v] * numRegions).T.flatten())


def smooth_trans(u0, u1, t, tc, teps):
    """
    Smooth transition from u0 to u1 in interval `tc - teps < t < tc + teps`.
    """
    t0 = tc - teps
    t1 = tc + teps
    return u0 if t <= t0 else u1 if t > t1 else \
        u0 + (u1 - u0) * (1 - np.cos(np.pi/(t1 - t0)*(t - t0))) * 0.5


class Ode:
    params_fixed = dict()
    params_prior = dict()
    """
    Solves model equations.
    params: dict()
        Parameters.
    t_span: `2-tuple of floats`
        Interval of integration.
    y0: `array_like`, (n_vars, n_regions)
        Initial state.
    t_eval: `array_like`, (nt)
        Times at which to store the computed solution.
    Returns:
    y: `array_like`, (n_vars, n_regions)
        Solution.
    """
    def solve(self, params, t_span, y0, t_eval):
        raise NotImplementedError()

    """
    Solves model equations converting initial state and solution from/to [S, I]
    S (susceptible), I (infected)

    params: dict()
        Parameters.
    t_span: `2-tuple of floats`
        Interval of integration.
    si0: `array_like`, (2, n_regions)
        Initial state.
    t_eval: `array_like`, (nt)
        Times at which to store the computed solution.
    Returns:
    si: `array_like`, (2, n_regions, nt)
        Solution.
    """

    def solve_SI(self, params, t_span, y0, t_eval):
        raise NotImplementedError()


class Data:
    total_infected = None
    time = None
    population = None
    names = None
    commute_matrix = None


class Sir(Ode):
    params_fixed = {'R0': 1.002, 'gamma': 60.}
    params_prior = {'R0': (0.5, 3), 'gamma': (1, 100)}

    def solve(self, params, t_span, y0, t_eval):
        """
        y0: `array_like`, (2, n_regions)
        """
        def rhs(t, y_flat):
            S, I = y_flat.reshape(2, -1)
            N = params['N']
            gamma = params['gamma']
            beta = params['R0'] * gamma
            c1 = beta * S * I / N
            c2 = gamma * I
            dSdt = -c1
            dIdt = c1 - c2
            return np.array((dSdt, dIdt)).flatten()

        y0 = np.array(y0)
        n_vars = 2
        assert y0.shape[0] == n_vars
        n_regions = y0.shape[1]
        y0_flat = y0.flatten()
        sol = solve_ivp(rhs, t_span=t_span, y0=y0_flat, t_eval=t_eval)
        y_flat = sol.y
        y = y_flat.reshape(n_vars, n_regions, len(t_eval))
        return y

    def solve_SI(self, params, t_span, si0, t_eval):
        si = self.solve(params, t_span, si0, t_eval)
        return si


class Seir(Ode):
    params_fixed = {
        'R0': 1.75,
        'Z': 2,
        'D': 0.8,
        'tact': 24.,
        'kbeta': 0.5,
        'nu': 1,
    }
    params_prior = {
        'R0': (0.5, 4),
        'Z': (0.01, 10),
        'D': (0.01, 10),
        'tact': (0, 60),
        'kbeta': (0., 1.),
        'nu': (0.1, 10.),
    }

    def solve(self, params, t_span, y0, t_eval):
        def rhs(t, y_flat):
            S, E, I = y_flat.reshape(3, -1)
            N = params['N']
            C = params['C']
            Z = params['Z']
            D = params['D']
            beta = params['R0'] / D
            beta = smooth_trans(beta, beta * params['kbeta'], t,
                                params['tact'], 0)
            nu = params['nu']
            k1 = I + nu * np.dot(C + C.T, I / N)
            c1 = beta * S * k1 / N
            c2 = E / Z
            dSdt = -c1
            dEdt = c1 - c2
            dIdt = c2 - I / D
            return np.array((dSdt, dEdt, dIdt)).flatten()

        y0 = np.array(y0)
        n_vars = 3
        assert y0.shape[0] == n_vars
        n_regions = y0.shape[1]
        y0_flat = y0.flatten()
        sol = solve_ivp(rhs, t_span=t_span, y0=y0_flat, t_eval=t_eval)
        y_flat = sol.y
        y = y_flat.reshape(n_vars, n_regions, len(t_eval))
        return y

    def solve_SI(self, params, t_span, si0, t_eval):
        S0, I0 = np.array(si0)
        E0 = np.zeros_like(S0)
        S, E, I = self.solve(params, t_span, [S0, E0, I0], t_eval)
        return np.array((S, I))


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
        I0 = Itotal[:, 0]
        S0 = N - I0
        y0 = S0, I0

        self.n_regions = Itotal.shape[0]
        self.region_names = data.name

        Idaily = np.diff(Itotal, axis=1)  # shape (n_regions, nt-1)
        Idaily = Idaily.T  # shape (nt-1, n_regions)

        self.data['Model']['x-data'] = repeat(t[1:], self.n_regions)
        self.data['Model']['y-data'] = Idaily.flatten()

        self.data['Model']['Initial Condition'] = y0
        self.data['Model']['Population Size'] = N

        if data.commute_matrix is None:
            C = np.zeros([self.n_regions] * 2)
        else:
            C = np.array(data.commute_matrix)
            assert list(C.shape) == [self.n_regions] * 2
        self.data['Model']['Commute Matrix'] = C

        T = np.ceil(t[-1])
        self.data['Propagation']['x-data'] = repeat(
            np.linspace(1, T, int(T + 1)), self.n_regions)

        save_file(self.data, self.saveInfo['inference data'],
                  'Data for Inference', 'pickle')

    def set_variables_and_distributions(self):
        p = self.params_to_infer + ['[r]']
        for k, name in enumerate(p):
            self.e['Variables'][k]['Name'] = name
            self.e['Variables'][k]['Prior Distribution'] = 'Prior for ' + name

        self.nParameters = len(p)

        k = 0

        for name in self.params_to_infer:
            self.e['Distributions'][k]['Name'] = 'Prior for ' + name
            self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
            minmax = self.ode.params_prior[name]
            self.e['Distributions'][k]['Minimum'] = minmax[0]
            self.e['Distributions'][k]['Maximum'] = minmax[1]
            k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for [r]'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 0.01
        self.e['Distributions'][k]['Maximum'] = 10.
        k += 1

    def get_params(self, korali_p, N, C):
        """
        N: `array_like`, (n_regions)
        Population of regions.
        C: `array_like`, (n_regions,n_regions)
        Commute matrix.
        """
        params = dict()
        params['N'] = np.array(N)
        params['C'] = np.array(C)
        for name, value in self.ode.params_fixed.items():
            params[name] = value
        for i, name in enumerate(self.params_to_infer):
            params[name] = korali_p[i]
        return params

    def computational_model(self, s):
        p = s['Parameters']
        t = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']
        C = self.data['Model']['Commute Matrix']

        tt1 = [min(t) - 1] + list(t[0::self.n_regions])
        assert min(tt1) >= 0

        params = self.get_params(p, N, C)
        S, _ = self.ode.solve_SI(params, [0, max(t)], y0, tt1)
        # S: shape (n_regions, nt)
        Idaily = -np.diff(S, axis=1)  # shape (n_regions, nt-1)
        Idaily = Idaily.T  # shape (nt-1, n_regions)
        Idaily = list(Idaily.flatten())

        s['Reference Evaluations'] = Idaily
        s['Dispersion'] = len(Idaily) * [p[-1]]

    def computational_model_propagate(self, s):
        p = s['Parameters']
        t = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']
        C = self.data['Model']['Commute Matrix']

        t1 = list(t[0::self.n_regions])
        assert min(t1) >= 0

        params = self.get_params(p, N, C)
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

        js['Dispersion'] = len(t1) * [p[-1]]
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
        self.propagatedVariables[varName] = np.zeros((Ns,Nt))
        for k in range(Ns):
            self.propagatedVariables[varName][k] = np.asarray(js['Dispersion'])


def moving_average(x, w):
    '''
    x: `numpy.ndarray`, (N)
    w: int
      Window half-width.
    Returns:
    xa: `numpy.ndarray`, (N)
      Array `x` averaged over window [-w,w].
    '''
    s = np.zeros_like(x)
    q = np.zeros_like(x)
    for i in range(len(x)):
        for j in range(max(0, i - w), min(i + w + 1, len(x))):
            if np.isfinite(x[j]):
                s[i] += x[j]
                q[i] += 1
    xa = s / q
    return xa


def fill_nans_nearest(x):
    '''
    x: `numpy.ndarray`, (N)
    Returns:
    `numpy.ndarray`, (N)
      Array `x` where each NaN value is replaced by a finite value from:
      - nearest index to the left
      - if not found, nearest index to the right
    '''
    def left(x, i):
        for v in reversed(x[:i + 1]):
            if np.isfinite(v):
                return v
        return x[i]

    def right(x, i):
        for v in x[i:]:
            if np.isfinite(v):
                return v
        return x[i]

    x = np.array(x)
    for i in range(len(x)):
        x[i] = left(x, i)
        x[i] = right(x, i)
    return x


def fill_nans_interp(t, x):
    from scipy.interpolate import interp1d
    x = np.copy(x)
    nans = np.isnan(x)

    #x[nans]= np.interp(t[nans], t[~nans], x[~nans], left=np.nan, right=np.nan)

    f = interp1d(t[~nans], x[~nans], fill_value='extrapolate')
    for i in range(len(x)):
        if np.isnan(x[i]):
            x[i] = f(t[i])
    return x


def get_data_switzerland() -> Data:
    data = Data()
    from epidemics.data.combined import RegionalData
    regionalData = RegionalData('switzerland')
    data.total_infected = regionalData.infected

    data.total_infected = moving_average(data.total_infected, 0)
    data.time = regionalData.time
    data.total_infected = data.total_infected
    data.time = data.time
    data.population = regionalData.populationSize
    return data


def get_data_switzerland_cantons(keys=None) -> Data:
    """
    keys: `array_like` or None
    List of canton keys to select (e.g ['ZH', 'TI']), select all if None.
    """
    key_to_total_infected = swiss_cantons.fetch_openzh_covid_data()
    key_to_population = swiss_cantons.CANTON_POPULATION

    if keys is None:
        keys = sorted(list(key_to_population))
    n_regions = len(keys)

    data = Data()
    if n_regions == 0:
        return data

    nt = len(key_to_total_infected[keys[0]])
    data.time = np.arange(nt)
    data.total_infected = np.empty((n_regions, nt))
    data.population = np.empty((n_regions))

    for i, k in enumerate(keys):
        Itotal = key_to_total_infected[k]
        #Itotal = fill_nans_nearest(Itotal)
        Itotal[0] = 1
        Itotal = fill_nans_interp(data.time, Itotal)
        Itotal = moving_average(Itotal, 2)
        data.total_infected[i, :] = Itotal[:]
        data.population[i] = key_to_population[k]
    data.name = keys
    return data


def get_data_synthetic(ode=Seir()) -> Data:
    data = Data()
    n_regions = 3
    N = np.array([8e6] * n_regions) / (10**np.arange(n_regions))
    I0 = np.array([16] * n_regions)
    t = np.arange(0, 60)
    params = {'N': N}
    params.update(ode.params_fixed)
    S0 = N - I0
    S, I = ode.solve_SI(params, [0, max(t)], [S0, I0], t)
    data.total_infected = N[:, None] - S
    data.time = t
    data.population = N
    return data


def main():
    nSamples = 1000
    x = argparse.Namespace()
    x.dataFolder = "data/"
    x.nPropagation = 20
    x.percentages = [0.5]
    x.nThreads = 8

    #data = get_data_switzerland()
    data = get_data_switzerland_cantons(['ZH', 'TI', 'VD'])
    N = data.population
    data.commute_matrix = np.array([
        [0, 0, 0],  #
        [0, 0, 0],  #
        [0, 0, 0],  #
    ])
    #data = get_data_synthetic()

    #ode = Sir()
    #params_to_infer = ['R0', 'gamma']

    ode = Seir()
    params_to_infer = ['R0', 'Z', 'D']

    a = Model(data, ode, params_to_infer, **vars(x))

    #a.sample(nSamples)
    #a.propagate()
    a.evaluate([1.8, 1., 2., 2.5])

    a.save()


if __name__ == "__main__":
    main()
