#!/usr/bin/env python3
import numpy as np

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data.swiss_cantons import CANTON_KEYS_ALPHABETICAL, CANTON_POPULATION
from epidemics.cantons.py.model import \
        get_canton_model_data, get_canton_reference_data, \
        get_municipality_model_data, ModelData

import libsolver

data = get_canton_model_data(include_foreign=False)
def example_run_seiin(num_days: int, inputs):
    """Runs the SEIIN model for some set of parameters and some initial conditions."""

    # Parameters.
    # Li2020, Table 1
    params = libsolver.solvers.seiin.Parameters(
            beta=inputs[0], mu=inputs[1], alpha=inputs[2], Z=inputs[3], D=inputs[4], theta=inputs[5])

    # Initial state.
    N0 = list(data.region_population)
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    IR0[data.key_to_index['TI']] = 1  # Ticino.
    IU0[data.key_to_index['TI']] = 100# We assume 100 unreported cases as well

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libsolver.solvers.seiin.Solver(data.to_cpp(), verbose=False)
    return solver.solve(params, y0, num_days)


            #All_results[isim][int(day/c)][:] = results[day].Iu()
            #reported   [isim][int(day/c)][:] = results[day].E()
