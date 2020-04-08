#!/usr/bin/env python3

"""Solve the ODE and plot."""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from data import CANTON_POPULATION, get_symmetric_Mij, fetch_canton_data
from plot import Renderer
from misc import Values, flatten

BUILD_DIR = os.path.join(os.path.dirname(__file__), '..', 'build')
if os.path.exists(BUILD_DIR):
    sys.path.append(BUILD_DIR)

import libsolver

CANTON_TO_INDEX, REFDATA = fetch_canton_data()
NUM_CANTONS = len(CANTON_TO_INDEX)

def example_run(num_days):
    """Runs the SEIIR model for some set of parameters and some initial conditions."""
    # Data.
    Mij = get_symmetric_Mij(CANTON_TO_INDEX)

    # Parameters.
    # Li2020, Table 1
    params = libsolver.Parameters (beta=1.12, mu=0., alpha=1., Z=3.69, D=3.47, theta=1.36)

    # Initial state.
    N0 = [CANTON_POPULATION[canton] for canton in CANTON_TO_INDEX]
    S0 = [CANTON_POPULATION[canton] for canton in CANTON_TO_INDEX]
    E0 = [0] * NUM_CANTONS
    IR0 = [0] * NUM_CANTONS
    IU0 = [0] * NUM_CANTONS
    IR0[CANTON_TO_INDEX['TI']] = 1  # Ticino.
    IU0[CANTON_TO_INDEX['TI']] = 0
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libsolver.Solver(flatten(Mij))
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
        for i, c in enumerate(rend.get_codes()):
            Ir = state.Ir(CANTON_TO_INDEX[c])
            print("{:02d} {} {:.1f}".format(i, c, Ir))
            values[c] = Ir / Ir_max * 2
            texts[c] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    rend = Renderer(frame_callback)
    rend.save_image()
    rend.save_movie(frames=len(results))


def plot_timeseries(results, cantons, var='S', refdata=False):
    VAR_NmS = "N-S"
    VAR_S0mS = "S0-S"
    if not isinstance(cantons, list):
        cantons = [cantons]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    ax.set_title("   ".join(
        ["$N_{{{:}}}$={:.0f}K".format(c, CANTON_POPULATION[c]*1e-3) for c in cantons]))
    vardesc = {
            'Ir' : "Infected Reported",
            'Iu' : "Infected Unreported",
            VAR_NmS : "N - S(t)",
            VAR_S0mS : "S(0) - S(t)",
            }
    ax.set_xlabel("days")
    ax.set_ylabel(vardesc[var])
    for c in cantons:
        if var == VAR_NmS:
            u = np.array([state.S(CANTON_TO_INDEX[c]) for state in results])
            u = CANTON_POPULATION[c] - u
        elif var == VAR_S0mS:
            u = np.array([state.S(CANTON_TO_INDEX[c]) for state in results])
            u = u[0] - u
        else:
            u = [getattr(state, var)(CANTON_TO_INDEX[c]) for state in results]
        line, = ax.plot(u,label=c)
        if refdata:
            u = [d[CANTON_TO_INDEX[c]] for d in REFDATA]
            ax.plot(u, c=line.get_color(), ls='', marker='.')

    if refdata:
        ax.set_ylim([1, 10 * np.nanmax([d[CANTON_TO_INDEX[c]] for c in cantons for d in REFDATA])])
    ax.set_yscale('log')
    ax.legend()
    varname = {
            'Ir' : "Ir",
            'Iu' : "Iu",
            VAR_NmS : "NmS",
            VAR_S0mS : "S0mS",
            }
    fig.tight_layout()
    fig.savefig("{:}_{:}.pdf".format(varname[var], "_".join(cantons)))


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('type', type=str, choices=('video', 'timeseries'), help="Plot type.")
    parser.add_argument('days', type=int, default=50, help="Number of days to evaluate.")
    args = parser.parse_args(argv)

    results = example_run(args.days)

    if args.type == 'video':
        plot_ode_results(results)
    elif args.type == 'timeseries':
        cc = ['TI', 'ZH', 'AG']
        #plot_timeseries(results, cc, var='Ir', refdata=True)
        plot_timeseries(results, cc, var="N-S", refdata=True)
        #plot_timeseries(results, cc, var="S0-S", refdata=True)


if __name__ == '__main__':
    main(sys.argv[1:])
