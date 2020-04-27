#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from epidemics.tools.tools import load_model

from pathlib import Path
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import time

# required for load_model
from main import Model, Ode, Sir, Seir, SeirCpp

from epidemics.tools.tools import prepare_folder


def compute_plot_intervals(model, varName, ns, ax, ylabel, cummulate=-1):
    xdata = model.data['Propagation']['x-data'][::model.n_regions]
    Ns = model.propagatedVariables[varName].shape[0]
    Nt = model.propagatedVariables[varName].shape[1]

    samples = np.zeros((Ns * ns, Nt))

    print(
        f"[Epidemics] Sampling from {model.likelihoodModel} for '{varName}' variable... ",
        end='',
        flush=True)

    start = time.process_time()

    if model.likelihoodModel == 'Normal':
        for k in range(Nt):
            m = model.propagatedVariables[varName][:, k]
            r = model.propagatedVariables['Standard Deviation'][:, k]
            x = [np.random.normal(m, r) for _ in range(ns)]
            samples[:, k] = np.asarray(x).flatten()

    elif model.likelihoodModel == 'Positive Normal':
        for k in range(Nt):
            m = model.propagatedVariables[varName][:, k]
            s = model.propagatedVariables['Standard Deviation'][:, k]
            t = get_truncated_normal(m, s, 0, np.Inf)
            x = [t.rvs() for _ in range(ns)]
            samples[:, k] = np.asarray(x).flatten()

    elif model.likelihoodModel == 'Negative Binomial':
        for k in range(Nt):
            m = model.propagatedVariables[varName][:, k]
            r = model.propagatedVariables['Dispersion'][:, k]
            p = p = m / (m + r)
            x = [np.random.negative_binomial(r, 1 - p) for _ in range(ns)]
            samples[:, k] = np.asarray(x).flatten()

    else:
        sys.exit(
            "\n[Epidemics] Likelihood not found in compute_plot_intervals.\n")

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

    for p in np.sort(model.percentages)[::-1]:
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
    if (model.logPlot): ax.set_yscale('log')

    plt.draw()

    return samples


def plot_intervals(model, region=0):
    if region >= model.n_regions:
        print("skpping region={:}".format(region))
        return
    print('[Epidemics] Compute and Plot credible intervals.')
    fig = plt.figure(figsize=(12, 8))
    if model.region_names:
        fig.suptitle(model.region_names[region])

    xdata = model.data['Model']['x-data'][region::model.n_regions]
    ydata = model.data['Model']['y-data'][region::model.n_regions]

    ax = fig.subplots(2)
    ax[0].plot(xdata,
               ydata,
               'o',
               lw=2,
               label='Daily Infected(data)',
               color='black')

    var = 'Daily Incidence {:}'.format(region)

    compute_plot_intervals(model, var, model.nPropagation, ax[0],
                           'Daily Incidence')

    z = np.cumsum(ydata)
    ax[1].plot(xdata,
               z,
               'o',
               lw=2,
               label='Cummulative Infected (data)',
               color='black')

    compute_plot_intervals(model,
                           var,
                           model.nPropagation,
                           ax[1],
                           'Cummulative number of infected',
                           cummulate=1)

    ax[-1].set_xlabel('time in days')

    name = "prediction{:}.png".format(str(region) if region else "")
    fig.savefig(name)
    plt.close(fig)

def plot_all_regions(model, names=None):
    print('[Epidemics] Plot data of all regions.')

    if names is None:
        names = model.region_names

    fig = plt.figure(figsize=(10, 10))
    axes = fig.subplots(6,5)
    axes = axes.flatten()

    iax = 0
    imax = min(len(axes), len(names))
    for ax,name in zip(axes[:imax], names[:imax]):
        region = model.region_names.index(name)
        ax.set_title(names[region], loc='left')

        # data
        x = model.data['Model']['x-data'][region::model.n_regions]
        y = model.data['Model']['y-data'][region::model.n_regions]
        ycum = np.cumsum(y)
        ax.plot(x, ycum, 'o', lw=1, color='black')
        ax.set_ylim(0, 60)
        ax.set_xticks([0, 30, 60])
        ax.set_ylim(0, ycum.max() * 1.5)
        ax.set_yticks([0, int(ycum.max())])

        # inference
        var = "Daily Incidence {:}".format(region)
        Ns = model.propagatedVariables[var].shape[0]
        Nt = model.propagatedVariables[var].shape[1]
        ns = 1
        samples = np.zeros((Ns * ns, Nt))
        if model.likelihoodModel == 'Negative Binomial':
            for k in range(Nt):
                m = model.propagatedVariables[var][:, k]
                r = model.propagatedVariables['Dispersion'][:, k]
                p = p = m / (m + r)
                y = [np.random.negative_binomial(r, 1 - p) for _ in range(ns)]
                samples[:, k] = np.asarray(y).flatten()
        samples = np.cumsum(samples, axis=1)
        mean = np.zeros((Nt, 1))
        median = np.zeros((Nt, 1))
        for k in range(Nt):
            median[k] = np.quantile(samples[:, k], 0.5)
            mean[k] = np.mean(samples[:, k])
        # one sample
        y = model.propagatedVariables[var][0, :]
        ycum = np.cumsum(y)
        x = model.data['Propagation']['x-data'][::model.n_regions]
        ax.plot(x, mean, lw=1, color='red')
        ax.plot(x, ycum, lw=1, color='blue')

    for ax in axes[imax:]:
        ax.set_axis_off()

    name = "all.png".format(str(region) if region else "")
    fig.tight_layout()
    fig.savefig(name)
    plt.close(fig)

def main():
    dataFolder = Path("data")
    f = dataFolder / 'cantons' / 'state.pickle'

    assert os.path.isfile(f)

    model = load_model(f)

    #for region in range(model.n_regions):
    #    plot_intervals(model, region=region)

    # plot data
    plot_all_regions(model)


if __name__ == "__main__":
    main()
