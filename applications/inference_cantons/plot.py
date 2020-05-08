#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))


def printlog(msg):
    out = sys.stderr
    out.write(str(msg) + "\n")
    out.flush()


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

"""
Always returns `int` (while `round(np.float32(1))` returns `np.float32`)
"""
def intround(a):
    return round(float(a))


def compute_plot_intervals(model, varName, ns, ax, ylabel, cummulate=-1):
    xdata = model.data['Propagation']['x-data'][::model.n_regions]
    Ns = model.propagatedVariables[varName].shape[0]
    Nt = model.propagatedVariables[varName].shape[1]

    samples = np.zeros((Ns * ns, Nt))

    printlog(
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
    printlog(f" elapsed {elapsed:.2f} sec")

    printlog(f"[Epidemics] Computing quantiles... ")

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
        printlog("skpping region={:}".format(region))
        return
    printlog('[Epidemics] Compute and Plot credible intervals.')
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


def plot_tiles(model, keys=None, daily=False):
    printlog('[Epidemics] Plot tiles with all regions.')

    cumulative = not daily

    if keys is None:
        keys = model.region_names

    keys = sorted(keys,
                  key=lambda key: np.cumsum(model.data['Model']['y-data'][
                      model.region_names.index(key)::model.n_regions]).max())

    #fig = plt.figure(figsize=(10, 10))
    #axes = fig.subplots(6, 5)
    fig = plt.figure(figsize=(20, 3.5))
    axes = fig.subplots(2, 13)
    axes = axes.flatten()

    iax = 0
    imax = min(len(axes), len(keys))
    for ax, key in zip(axes[:imax], keys[:imax]):
        region = model.region_names.index(key)
        ax.set_title(key, loc='left')

        # data
        x = model.data['Model']['x-data'][region::model.n_regions]
        y = model.data['Model']['y-data'][region::model.n_regions]
        ycum = np.cumsum(y)
        x_days = np.linspace(1, intround(max(x)), intround(max(x)))
        ycum_days = np.interp(x_days, x, ycum)
        y_days = np.diff(ycum_days, prepend=0)
        xplot = x_days
        yplot = np.cumsum(y_days) if cumulative else y_days
        ax.scatter(xplot, yplot, s=1, color='black')
        ax.set_xticks([0, max(xplot) // 2, max(xplot)])
        ax.set_ylim(0, yplot.max() * 2)
        ax.set_yticks([0, int(yplot.max())])

        # mean and median prediction
        x = model.data['Propagation']['x-data'][::model.n_regions]
        var = "Daily Incidence {:}".format(region)
        Ns = model.propagatedVariables[var].shape[0]
        Nt = model.propagatedVariables[var].shape[1]
        ns = 5
        samples = np.zeros((Ns * ns, Nt))
        if model.likelihoodModel == 'Negative Binomial':
            for k in range(Nt):
                m = model.propagatedVariables[var][:, k]
                r = model.propagatedVariables['Dispersion'][:, k]
                p = p = m / (m + r)
                y = [np.random.negative_binomial(r, 1 - p) for _ in range(ns)]
                samples[:, k] = np.asarray(y).flatten()
        if cumulative:
            samples = np.cumsum(samples, axis=1)
        mean = np.zeros((Nt, 1))
        median = np.zeros((Nt, 1))
        for k in range(Nt):
            median[k] = np.quantile(samples[:, k], 0.5)
            mean[k] = np.mean(samples[:, k])
        ax.plot(x, mean, lw=1, color='red')

        # one sample
        y = model.propagatedVariables[var][0, :]
        yplot = np.cumsum(y) if cumulative else y
        #ax.plot(x, ycum, lw=1, color='blue')

    for ax in axes[imax:]:
        ax.set_axis_off()

    fig.tight_layout()
    fn = "tiles_daily.png" if daily else "tiles.png"
    printlog(fn)
    fig.savefig(fn)
    plt.close(fig)


def get_data(model, keys):
    """
    model: `json`
        Korali model.
    keys: `array_likt`
        List of region keys to extract.
    Returns reference data for infections saved in `model`.
    t: `numpy.ndarray`, (nt)
        Days.
    I: `dict(key : numpy.ndarray)`, (nt)
        Daily infections.
    Itotal: `dict(key : numpy.ndarray)`, (nt)
        Total number of infections.
    """
    t = None
    I = dict()
    Itotal = dict()
    for key in keys:
        region = model.region_names.index(key)
        t = model.data['Model']['x-data'][region::model.n_regions]
        y = model.data['Model']['y-data'][region::model.n_regions]
        y = np.array(y).flatten()
        I[key] = y
        Itotal[key] = np.cumsum(y)
    return t, I, Itotal


def get_fit(model, keys, n_samples=5):
    """
    model: `json`
        Korali model.
    keys: `array_likt`
        List of region keys to extract.
    Returns reference data for infections saved in `model`.
    t: `numpy.ndarray`, (nt)
        Days.
    I: `dict(key : numpy.ndarray)`, (nt)
        Mean fit for daily infections.
    Itotal: `dict(key : numpy.ndarray)`, (nt)
        Mean fit for total number of infections.
    """
    t = None
    I = dict()
    Itotal = dict()
    for key in keys:
        region = model.region_names.index(key)
        t_data = model.data['Model']['x-data'][region::model.n_regions]
        t_data = np.array(t_data).flatten()

        var = "Daily Incidence {:}".format(region)
        Ns = model.propagatedVariables[var].shape[0]
        Nt = model.propagatedVariables[var].shape[1]
        t = np.linspace(0, max(t_data), Nt)
        ns = n_samples
        samples = np.zeros((Ns * ns, Nt))
        if model.likelihoodModel == 'Negative Binomial':
            for it in range(Nt):
                m = model.propagatedVariables[var][:, it]
                r = model.propagatedVariables['Dispersion'][:, it]
                p = p = m / (m + r)
                y = [np.random.negative_binomial(r, 1 - p) for _ in range(ns)]
                samples[:, it] = np.asarray(y).flatten()
        else:
            raise NotImplementedError()

        samples_cum = np.cumsum(samples, axis=1)
        mean = np.zeros(Nt)
        mean_cum = np.zeros(Nt)
        for it in range(Nt):
            mean[it] = np.mean(samples[:, it])
            mean_cum[it] = np.mean(samples_cum[:, it])

        I[key] = mean
        Itotal[key] = mean_cum
    return t, I, Itotal


class MapType:
    data = "data"
    fit = "fit"
    error = "error"


def plot_map(model,
             map_type=MapType.fit,
             image=True,
             movie=True,
             image_day=-1):
    keys = model.region_names
    nt = None

    if map_type == MapType.data:
        t, _, Itotal = get_data(model, keys)
        nt = len(t)
    elif map_type == MapType.fit:
        t, _, Itotal = get_fit(model, keys)
        nt = len(t)
    elif map_type == MapType.error:
        t_data, _, Itotal_data = get_data(model, keys)
        t_fit, _, Itotal_fit = get_fit(model, keys)
        nt = intround(max(max(t_data), max(t_fit)))
        t = np.linspace(1, nt, nt)

        Itotal_error = dict()
        for key in keys:
            Id = np.interp(t, t_data, Itotal_data[key])
            If = np.interp(t, t_fit, Itotal_fit[key])
            Itotal_error[key] = (If - Id) / np.maximum(1, Id)
    else:
        raise NotImplementedError()

    from epidemics.cantons.py.plot import Renderer

    def frame_callback(rend):
        values = dict()
        texts = dict()
        colors = dict()
        for c in rend.get_codes():
            i = (nt - 1) * rend.get_frame() // rend.get_max_frame()
            if map_type in [MapType.data, MapType.fit]:
                v = Itotal[c][i]
                values[c] = v * 1e-3
                texts[c] = "{:}".format(int(v))
            if map_type in [MapType.error]:
                cmap = plt.get_cmap("coolwarm")
                v = Itotal_error[c][i]
                vclip = np.clip(v, -1, 1)
                values[c] = abs(vclip) * 0.5
                texts[c] = "{:+d}%".format(int(v * 100))
                colors[c] = cmap(vclip * 0.5 + 0.5)
        rend.set_values(values)
        rend.set_texts(texts)
        rend.set_colors(colors)

    from epidemics.cantons.py.model import get_canton_model_data
    rend = Renderer(frame_callback,
                    data=get_canton_model_data(),
                    resolution=(1920, 1080),
                    draw_Mij=False,
                    draw_Cij=True)

    fnbase = map_type
    frames = intround(max(t)) + 1
    if image:
        frame = min(image_day, frames - 1) if image_day != -1 else frames - 1
        fn = "{:}_day{:}.png".format(fnbase, frame)
        rend.save_image(frame=image_day, frames=frames, filename=fn)
        printlog(fn)
    if movie:
        fn = fnbase + ".mp4"
        rend.save_movie(frames=intround(max(t)), filename=fn)
        printlog(fn)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tiles',
                        action='store_true',
                        help="Tiles with cumulative infections in all regions")
    parser.add_argument('--tiles_daily',
                        action='store_true',
                        help="Tiles with daily new infections in all regions")
    parser.add_argument('--map_data',
                        action='store_true',
                        help="Map with data")
    parser.add_argument('--map_fit', action='store_true', help="Map with fit")
    parser.add_argument('--map_error',
                        action='store_true',
                        help="Map with error")
    parser.add_argument('--image', action='store_true', help="Create image")
    parser.add_argument('--image_day',
                        type=int,
                        default=-1,
                        help="Select day for image")
    parser.add_argument('--movie', action='store_true', help="Create movie")
    parser.add_argument('--intervals',
                        action='store_true',
                        help="Plot credible intervals for all regions")
    args = parser.parse_args()

    dataFolder = Path("data")
    f = dataFolder / 'cantons' / 'state.pickle'

    model = load_model(f)

    if args.intervals:
        for region in range(model.n_regions):
            plot_intervals(model, region=region)

    if args.tiles:
        plot_tiles(model)

    if args.tiles_daily:
        plot_tiles(model, daily=True)

    if args.map_data:
        plot_map(model,
                 map_type=MapType.data,
                 image=args.image,
                 movie=args.movie,
                 image_day=args.image_day)

    if args.map_fit:
        plot_map(model,
                 map_type=MapType.fit,
                 image=args.image,
                 movie=args.movie,
                 image_day=args.image_day)

    if args.map_error:
        plot_map(model,
                 map_type=MapType.error,
                 image=args.image,
                 movie=args.movie,
                 image_day=args.image_day)


if __name__ == "__main__":
    main()
