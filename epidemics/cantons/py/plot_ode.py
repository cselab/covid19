#!/usr/bin/env python3

"""Solve the ODE and plot the results."""

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data.swiss_cantons import CANTON_KEYS_ALPHABETICAL, CANTON_POPULATION
from epidemics.cantons.py.model import get_canton_model_data, get_canton_validation_data
from epidemics.cantons.py.plot import Renderer

import libsolver

KEY_TO_INDEX = {key: k for k, key in enumerate(CANTON_KEYS_ALPHABETICAL)}
KEY_TO_POPULATION = CANTON_POPULATION
VALIDATION_DATA = get_canton_validation_data()

def example_run(num_days, include_foreign=True):
    """Runs the SEIIR model for some set of parameters and some initial conditions."""
    # Start date is needed to compute properly the inflow of foreign infected people.
    model_data = get_canton_model_data(include_foreign=include_foreign)

    # Parameters.
    # Li2020, Table 1
    params = libsolver.Parameters (beta=1.12, mu=0., alpha=1., Z=3.69, D=3.47, theta=1.36)

    # Initial state.
    N0 = [KEY_TO_POPULATION[key] for key in KEY_TO_INDEX]
    E0 = [0] * len(KEY_TO_INDEX)
    IR0 = [0] * len(KEY_TO_INDEX)
    IU0 = [0] * len(KEY_TO_INDEX)
    IR0[KEY_TO_INDEX['TI']] = 1  # Ticino.
    IU0[KEY_TO_INDEX['TI']] = 0

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libsolver.Solver(model_data.to_cpp())
    results = solver.solve(params, y0, num_days)
    return results


def plot_ode_results(results):
    """Plot results from the ODE.

    Arguments:
        results: A list of State objects.
    """
    Ir_max = np.max([state.Ir() for state in results])

    def frame_callback(rend):
        t = rend.get_frame() * (len(results) - 1) // rend.get_max_frame()

        state = results[t]
        values = {}
        texts = {}
        for i, key in enumerate(rend.get_codes()):
            Ir = state.Ir(KEY_TO_INDEX[key])
            print("{:02d} {} {:.1f}".format(i, key, Ir))
            values[key] = Ir / Ir_max * 2
            texts[key] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    rend = Renderer(frame_callback)
    rend.save_image()
    rend.save_movie(frames=len(results))


def plot_timeseries(results, keys, var='S', refdata=False):
    VAR_NmS = "N-S"
    VAR_S0mS = "S0-S"
    if not isinstance(keys, list):
        keys = [keys]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    ax.set_title("   ".join(
        ["$N_{{{:}}}$={:.0f}K".format(key, KEY_TO_POPULATION[key]*1e-3) for key in keys]))
    vardesc = {
            'Ir' : "Infected Reported",
            'Iu' : "Infected Unreported",
            VAR_NmS : "N - S(t)",
            VAR_S0mS : "S(0) - S(t)",
            }
    ax.set_xlabel("days")
    ax.set_ylabel(vardesc[var])
    for key in keys:
        if var == VAR_NmS:
            u = np.array([state.S(KEY_TO_INDEX[key]) for state in results])
            u = KEY_TO_POPULATION[key] - u
        elif var == VAR_S0mS:
            u = np.array([state.S(KEY_TO_INDEX[key]) for state in results])
            u = u[0] - u
        else:
            u = [getattr(state, var)(KEY_TO_INDEX[key]) for state in results]
        line, = ax.plot(u,label=key)
        if refdata:
            u = VALIDATION_DATA.cases_per_country[key]
            ax.plot(u, c=line.get_color(), ls='', marker='.')

    if refdata:
        ax.set_ylim([1, 10 * np.nanmax(VALIDATION_DATA.cases_per_country[key])])
    ax.set_yscale('log')
    ax.legend()
    varname = {
            'Ir' : "Ir",
            'Iu' : "Iu",
            VAR_NmS : "NmS",
            VAR_S0mS : "S0mS",
            }
    fig.tight_layout()
    fig.savefig("{:}_{:}.pdf".format(varname[var], "_".join(keys)))


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('type', type=str, choices=('video', 'timeseries'), help="Plot type.")
    parser.add_argument('days', type=int, default=50, help="Number of days to evaluate.")
    parser.add_argument('--no-foreign', action='store_true', help="Disable foreign commuters from the model.")
    args = parser.parse_args(argv)

    results = example_run(args.days, include_foreign=not args.no_foreign)

    if args.type == 'video':
        plot_ode_results(results)
    elif args.type == 'timeseries':
        keys = ['TI', 'ZH', 'AG']
        #plot_timeseries(results, keys, var='Ir', refdata=True)
        plot_timeseries(results, keys, var="N-S", refdata=True)
        #plot_timeseries(results, keys, var="S0-S", refdata=True)


if __name__ == '__main__':
    main(sys.argv[1:])
