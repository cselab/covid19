#!/usr/bin/env python3

"""Solve the ODE and plot the results."""

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data.swiss_cantons import CANTON_KEYS_ALPHABETICAL, CANTON_POPULATION
from epidemics.cantons.py.model import \
        get_canton_model_data, get_canton_reference_data, \
        get_municipality_model_data, ModelData
from epidemics.cantons.py.plot import Renderer

import libsolver

def example_run_seiin(data: ModelData, num_days: int):
    """Runs the SEIIN model for some set of parameters and some initial conditions."""

    # Parameters.
    # Li2020, Table 1
    params = libsolver.solvers.seiin.Parameters(
            beta=1.12, mu=0., alpha=1., Z=3.69, D=3.47, theta=1.36)

    # Initial state.
    N0 = list(data.region_population)
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    if 'TI' in data.key_to_index:
        IR0[data.key_to_index['TI']] = 1  # Ticino.
        IU0[data.key_to_index['TI']] = 0
    else:
        IR0[data.key_to_index['MUN-5192']] = 1  # Lugano.
        IU0[data.key_to_index['MUN-5192']] = 0

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libsolver.solvers.seiin.Solver(data.to_cpp(), verbose=True)
    return solver.solve(params, y0, num_days)


def example_run_seii_c(data, num_days):
    """Runs the SEII_C model for some set of parameters and some initial conditions."""
    # Parameters.
    params = libsolver.solvers.seii_c.Parameters(
            beta=3.0, nu=3.0, alpha=0.6, Z=3.69, D=3.47)

    # Initial state.
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    # if 'TI' in data.key_to_index:
    #     IU0[data.key_to_index['TI']] = 10  # Ticino.
    # else:
    #     IU0[data.key_to_index['MUN-5192']] = 10  # Lugano.

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(data.region_population, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0

    # Run the ODE solver.
    solver = libsolver.solvers.seii_c.Solver(data.to_cpp(), verbose=True)
    return solver.solve(params, y0, num_days)


def plot_ode_results(data: ModelData, results):
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
        for i, key in enumerate(data.region_keys):
            Ir = state.Ir(data.key_to_index[key])
            Iu = state.Iu(data.key_to_index[key])
            print("{:02d} {} Ir={:4.1f} Iu={:4.1f}".format(i, key, Ir, Iu))
            values[key] = Ir / Ir_max * 2
            texts[key] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    rend = Renderer(frame_callback, data=data)
    rend.save_image()
    rend.save_movie(frames=len(results))


def plot_timeseries(data: ModelData, results, keys, var='S', ref_data=None):
    key_to_population = dict(zip(data.region_keys, data.region_population))

    VAR_NmS = "N-S"
    VAR_S0mS = "S0-S"
    if not isinstance(keys, list):
        keys = [keys]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    ax.set_title("   ".join(
        ["$N_{{{:}}}$={:.0f}K".format(key, key_to_population[key]*1e-3) for key in keys]))
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
            u = np.array([state.S(data.key_to_index[key]) for state in results])
            u = key_to_population[key] - u
        elif var == VAR_S0mS:
            u = np.array([state.S(data.key_to_index[key]) for state in results])
            u = u[0] - u
        else:
            u = [getattr(state, var)(data.key_to_index[key]) for state in results]
        line, = ax.plot(u,label=key)
        if ref_data:
            u = ref_data.cases_per_country[key]
            ax.plot(u, c=line.get_color(), ls='', marker='.')

    if ref_data:
        ax.set_ylim([1, 10 * np.nanmax(ref_data.cases_per_country[key])])
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
    parser.add_argument('--level', type=str, choices=('canton', 'municipality'), default='canton', help="Level of details.")
    parser.add_argument('--model', type=str, choices=('seiin', 'seii_c'), default='seiin', help="Model.")
    args = parser.parse_args(argv)

    if args.level == 'canton':
        model_data = get_canton_model_data(include_foreign=not args.no_foreign)
        ref_data = get_canton_reference_data()
    else:
        model_data = get_municipality_model_data()
        ref_data = None

    if args.model == 'seiin':
        results = example_run_seiin(model_data, args.days)
    else:
        model_data.Mij *= 0.0
        results = example_run_seii_c(model_data, args.days)

    if args.type == 'video':
        plot_ode_results(model_data, results)
    elif args.type == 'timeseries':
        keys = ['TI', 'ZH', 'AG']
        #plot_timeseries(model_data, results, keys, var='Ir', ref_data=ref_data)
        plot_timeseries(model_data, results, keys, var="N-S", ref_data=ref_data)
        #plot_timeseries(model_data, results, keys, var="S0-S", ref_data=ref_data)


if __name__ == '__main__':
    main(sys.argv[1:])
