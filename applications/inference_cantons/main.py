#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import copy
from epidemics.tools.tools import import_from
import argparse

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.tools.tools import prepare_folder, save_file
from epidemics.epidemics import EpidemicsBase
import time

def repeat(v, numRegions):
    return list(np.array([v] * numRegions).T.flatten())


def mul2(v, nr):
    if nr == 2:
        v = np.array(v).astype(float)
        v[::2] *= 0.5
        return list(v)
    return v


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
    populationSize = None

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
    params_fixed = {'R0': 1.75, 'Z': 2, 'D': 0.8, 'tact': 24., 'kbeta': 0.5}
    params_prior = {
        'R0': (0.5, 4),
        'Z': (0.01, 10),
        'D': (0.01, 10),
        'tact': (0, 60),
        'kbeta': (0., 1.)
    }

    def solve(self, params, t_span, y0, t_eval):
        def rhs(t, y_flat):
            S, E, I = y_flat.reshape(3, -1)
            N = params['N']
            Z = params['Z']
            D = params['D']
            beta = params['R0'] / D
            beta = smooth_trans(beta, beta * params['kbeta'], t,
                                params['tact'], 0)
            c1 = beta * S * I / N
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
        Itotal = np.array(data.total_infected) # shape (n_regions, nt)
        t = np.array(data.time)
        N = np.array(data.populationSize)
        I0 = Itotal[:,0]
        S0 = N - I0
        y0 = S0, I0

        self.n_regions = Itotal.shape[0]

        Idaily = np.diff(Itotal, axis=1)  # shape (n_regions, nt-1)
        Idaily = Idaily.T # shape (nt-1, n_regions)

        self.data['Model']['x-data'] = repeat(t[1:], self.n_regions)
        self.data['Model']['y-data'] = Idaily.flatten()

        self.data['Model']['Initial Condition'] = y0
        self.data['Model']['Population Size'] = N

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

    def get_params(self, korali_p, N):
        """
        N: `array_like`, (n_regions)
        Population of regions.
        """
        params = {'N': np.array(N)}
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

        tt1 = [min(t) - 1] + list(t[0::self.n_regions])
        assert min(tt1) >= 0

        params = self.get_params(p, N)
        S,_ = self.ode.solve_SI(params, [0, max(t)], y0, tt1)
        # S: shape (n_regions, nt)
        Idaily = -np.diff(S, axis=1)  # shape (n_regions, nt-1)
        Idaily = Idaily.T   # shape (nt-1, n_regions)
        Idaily = list(Idaily.flatten())

        s['Reference Evaluations'] = Idaily
        s['Dispersion'] = len(Idaily) * [p[-1]]

    def computational_model_propagate(self, s):
        p = s['Parameters']
        t = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']
        N = self.data['Model']['Population Size']

        t1 = list(t[0::self.n_regions])
        assert min(t1) >= 0

        params = self.get_params(p, N)
        S,_ = self.ode.solve_SI(params, [0, max(t)], y0, t1)
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

    def compute_plot_intervals(self, varName, ns, ax, ylabel, cummulate=-1):
        xdata = self.data['Propagation']['x-data'][::self.n_regions]
        Ns = self.propagatedVariables[varName].shape[0]
        Nt = self.propagatedVariables[varName].shape[1]

        samples = np.zeros((Ns * ns, Nt))

        print(
            f"[Epidemics] Sampling from {self.likelihoodModel} for '{varName}' variable... ",
            end='',
            flush=True)

        start = time.process_time()

        if self.likelihoodModel == 'Normal':
            for k in range(Nt):
                m = self.propagatedVariables[varName][:, k]
                r = self.propagatedVariables['Standard Deviation'][:, k]
                x = [np.random.normal(m, r) for _ in range(ns)]
                samples[:, k] = np.asarray(x).flatten()

        elif self.likelihoodModel == 'Positive Normal':
            for k in range(Nt):
                m = self.propagatedVariables[varName][:, k]
                s = self.propagatedVariables['Standard Deviation'][:, k]
                t = get_truncated_normal(m, s, 0, np.Inf)
                x = [t.rvs() for _ in range(ns)]
                samples[:, k] = np.asarray(x).flatten()

        elif self.likelihoodModel == 'Negative Binomial':
            for k in range(Nt):
                m = self.propagatedVariables[varName][:, k]
                r = self.propagatedVariables['Dispersion'][:, k]
                p = p = m / (m + r)
                x = [np.random.negative_binomial(r, 1 - p) for _ in range(ns)]
                samples[:, k] = np.asarray(x).flatten()

        else:
            sys.exit(
                "\n[Epidemics] Likelihood not found in compute_plot_intervals.\n"
            )

        if cummulate > 0:
            samples = np.cumsum(samples, axis=cummulate)

        elapsed = time.process_time() - start
        print(f" elapsed {elapsed:.2f} sec")

        print(f"[Epidemics] Computing quantiles... ")

        mean = np.zeros((Nt, 1))
        median = np.zeros((Nt, 1))
        for k in range(Nt):
            median[k] = np.quantile(samples[:, k], 0.5)
            mean[k] = np.mean(samples[:, k])

        for p in np.sort(self.percentages)[::-1]:
            q1 = np.zeros((Nt, ))
            q2 = np.zeros((Nt, ))
            for k in range(Nt):
                q1[k] = np.quantile(samples[:, k], 0.5 - p / 2)
                q2[k] = np.quantile(samples[:, k], 0.5 + p / 2)
            ax.fill_between(xdata,
                            q1,
                            q2,
                            alpha=0.5,
                            label=f' {100*p:.1f}% credible interval')

        ax.plot(xdata, mean, '-', lw=2, label='Mean', color='black')
        ax.plot(xdata, median, '--', lw=2, label='Median', color='blue')

        ax.legend(loc='upper left')
        ax.set_ylabel(ylabel)
        #ax.set_xticks( range( np.ceil( max( xdata )+1 ).astype(int) ) )
        ax.grid()
        if (self.logPlot): ax.set_yscale('log')

        plt.draw()

        return samples

    def plot_intervals(self, region=0):
        if region >= self.n_regions:
            print("skpping region={:}".format(region))
            return
        print('[Epidemics] Compute and Plot credible intervals.')
        fig = plt.figure(figsize=(12, 8))
        #fig.suptitle(self.modelDescription)

        xdata = self.data['Model']['x-data'][region::self.n_regions]
        ydata = self.data['Model']['y-data'][region::self.n_regions]

        ax = fig.subplots(2)
        ax[0].plot(xdata,
                   ydata,
                   'o',
                   lw=2,
                   label='Daily Infected(data)',
                   color='black')

        var = 'Daily Incidence {:}'.format(region)

        self.compute_plot_intervals(var, self.nPropagation, ax[0],
                                    'Daily Incidence')

        z = np.cumsum(ydata)
        ax[1].plot(xdata,
                   z,
                   'o',
                   lw=2,
                   label='Cummulative Infected (data)',
                   color='black')

        self.compute_plot_intervals(var,
                                    self.nPropagation,
                                    ax[1],
                                    'Cummulative number of infected',
                                    cummulate=1)

        ax[-1].set_xlabel('time in days')

        name = "prediction{:}.png".format(str(region) if region else "")
        f = os.path.join(self.saveInfo['figures'], name)
        prepare_folder(os.path.dirname(f), clean=False)
        fig.savefig(f)

        plt.show()

        plt.close(fig)


def get_data_switzerland() -> Data:
    data = Data()
    from epidemics.data.combined import RegionalData
    regionalData = RegionalData('switzerland')
    data.total_infected = regionalData.infected

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
        return s / q

    data.total_infected = moving_average(data.total_infected, 0)
    data.time = regionalData.time
    sel = len(data.total_infected)
    #sel = 30
    data.total_infected = data.total_infected[:sel]
    data.time = data.time[:sel]
    data.populationSize = regionalData.populationSize
    return data


def get_data_synthetic(ode=Seir()) -> Data:
    data = Data()
    n_regions = 3
    N = np.array([8e6] * n_regions) / (10 ** np.arange(n_regions))
    I0 = np.array([16] * n_regions)
    t = np.arange(1, 60)
    params = {'N':N}
    params.update(ode.params_fixed)
    S0 = N - I0
    S,I = ode.solve_SI(params, [0, max(t)], [S0, I0], t)
    data.total_infected = N[:,None] - S
    data.time = t
    data.populationSize = N
    return data


def main():
    nSamples = 1000
    x = argparse.Namespace()
    x.dataFolder = "data/"
    x.nPropagation = 20
    x.percentages = [0.5]
    x.nThreads = 8

    #data = get_data_switzerland()
    data = get_data_synthetic()

    #ode = Sir()
    #params_to_infer = ['R0', 'gamma']

    ode = Seir()
    params_to_infer = ['R0', 'Z', 'D']

    a = Model(data, ode, params_to_infer, **vars(x))

    a.sample(nSamples)
    a.propagate()
    a.save()
    #a.plot_intervals()


if __name__ == "__main__":
    main()
