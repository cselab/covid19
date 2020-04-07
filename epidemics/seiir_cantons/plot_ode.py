#!/usr/bin/env python3

"""Solve the ODE and plot."""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from data import CANTON_POPULATION, get_symmetric_Mij, fetch_canton_data
from plot import Renderer
from misc import Values, flatten, extract_values_from_state

BUILD_DIR = os.path.join(os.path.dirname(__file__), 'build')
if os.path.exists(BUILD_DIR):
    sys.path.append(BUILD_DIR)

import libseiir

CANTON_TO_INDEX, _ = fetch_canton_data()
NUM_CANTONS = len(CANTON_TO_INDEX)

def example_run(num_days):
    """Runs the SEIIR model for some set of parameters and some initial conditions."""
    # Data.
    Mij = get_symmetric_Mij(CANTON_TO_INDEX)

    # Parameters.
    params = libseiir.Parameters(beta=1.1, mu=0., alpha=1., Z=4., D=4., theta=1.)

    # Initial state.
    N0 = [CANTON_POPULATION[canton] for canton in CANTON_TO_INDEX]
    S0 = N0
    E0 = [0] * NUM_CANTONS
    IR0 = [0] * NUM_CANTONS
    IU0 = [0] * NUM_CANTONS
    IR0[CANTON_TO_INDEX['TI']] = 1  # Ticino.
    IU0[CANTON_TO_INDEX['TI']] = 0
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libseiir.MultiSEIIR(flatten(Mij))
    results = solver.solve(params, y0, num_days)
    return results


def plot_ode_results(results):
    """Plot results from the ODE.

    `results` is a list of vector states.
    A vector state is a concatenated list of 5 x 26 elements:
        [S..., E..., Ir..., Iu..., N...]
    where
        S... represents the S value for 26 cantons, in the order specified by CANTON_TO_INDEX.
    """

    def frame_callback(rend):
        t = rend.get_frame() * (len(results) - 1) // rend.get_max_frame()
        state = results[t]
        Ir_per_canton = extract_values_from_state(state, NUM_CANTONS, Values.Ir)

        Ir_max = np.max([
            extract_values_from_state(state, NUM_CANTONS, Values.Ir)
            for state in results])

        values = dict()
        texts = dict()
        for i, c in enumerate(rend.get_codes()):
            Ir = Ir_per_canton[CANTON_TO_INDEX[c]]
            print("{:02d} {} {:.1f}".format(i, c, Ir))
            values[c] = Ir / Ir_max * 2
            texts[c] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    rend = Renderer(frame_callback)
    #rend.save_image()
    rend.save_movie(frames=len(results))

def plot_timeseries(results, cantons, var=Values.Ir):
    if not isinstance(cantons, list):
        cantons = [cantons]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    ax.set_title("   ".join(
        ["$N_{{{:}}}$={:.0f}K".format(c, CANTON_POPULATION[c]*1e-3) for c in cantons]))
    vardesc = {
            Values.Ir : "Infected Reported",
            Values.Iu : "Infected Unreported",
            }
    ax.set_xlabel("days")
    ax.set_ylabel(vardesc[var])
    for canton in cantons:
        u = [extract_values_from_state(state, NUM_CANTONS, var)
                [CANTON_TO_INDEX[canton]] for state in results]
        ax.plot(u,label=canton)
    ax.legend()
    varname = {
            Values.Ir : "Ir",
            Values.Iu : "Iu",
            }
    fig.tight_layout()
    fig.savefig("{:}_{:}.pdf".format(varname[var], "_".join(cantons)))

def main():
    num_days = 120
    results = example_run(num_days)
    #plot_ode_results(results)
    plot_timeseries(results, ['ZH', 'AG', 'TI'])


if __name__ == '__main__':
    main()
