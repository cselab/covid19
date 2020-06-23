#!/usr/bin/env python3

"""Solve the ODE for many different input parameters (needed for optimal sensor placement)"""

import numpy as np
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data.swiss_cantons import CANTON_KEYS_ALPHABETICAL, CANTON_POPULATION
from epidemics.cantons.py.model import \
        get_canton_model_data, get_canton_reference_data, \
        get_municipality_model_data, ModelData

import libepidemics

data = get_canton_model_data(include_foreign=False)
def example_run_seiin(num_days, inputs):
    num_days = int(num_days)

    cantons_ = [0,3,4,5,6,7,9,15,20,22,23,25]
    params = libepidemics.cantons.seiin_interventions.Parameters(
             beta =inputs[0],
             mu   =inputs[1],
             alpha=inputs[2],
             Z    =inputs[3],
             D    =inputs[4],
             theta=inputs[5],
             b1   =inputs[6],
             b2   =inputs[7],
             b3   =inputs[8],
             d1   =inputs[9],
             d2   =inputs[10],
             d3   =inputs[11])

    # Initial state.
    N0  = list(data.region_population)
    E0  = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions
    IR0[data.key_to_index['TI']] = 1  # Ticino.
    k = 0
    for c in cantons_:
        IU0[c] =  inputs[12+k]*N0[c]
        E0[c] = 3*inputs[12+k]*N0[c]
        k += 1

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libepidemics.cantons.seiin_interventions.Solver(data.to_cpp())
    y0 = libepidemics.cantons.seiin_interventions.State(y0)
    return solver.solve(params, y0, t_eval=range(1, num_days + 1))
