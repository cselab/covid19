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


    L = len(inputs)
    ic_cantons = len(cantons_)

    N0  = list(data.region_population)
    E0  = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions
    IR0[data.key_to_index['TI']] = 1  # Ticino.

    if L == 6 + ic_cantons: #case 1 

       params = libepidemics.cantons.seiin_interventions.Parameters(
                 beta  =inputs[0],
                 mu    =inputs[1],
                 alpha =inputs[2],
                 Z     =inputs[3],
                 D     =inputs[4],
                 theta =inputs[5],
                 b1    =inputs[0],
                 b2    =inputs[0],
                 b3    =inputs[0],
                 d1    =num_days+1,
                 d2    =num_days+1,
                 d3    =num_days+1,
                 theta1=inputs[5],
                 theta2=inputs[5],
                 theta3=inputs[5])

       k = 0
       for c in cantons_:
          IU0[c] =   inputs[6+k]
          E0 [c] = 3*inputs[6+k]
          k += 1

       S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
       y0 = S0 + E0 + IR0 + IU0 + N0

       # Run the ODE solver.
       solver = libepidemics.cantons.seiin_interventions.Solver(data.to_cpp())
       y0 = libepidemics.cantons.seiin_interventions.State(y0)
       return solver.solve(params, y0, t_eval=range(1, num_days + 1))


    if L == 8 + ic_cantons: #case 2 

       params = libepidemics.cantons.seiin_interventions.Parameters(
                 beta  =inputs[0],
                 mu    =inputs[1],
                 alpha =inputs[2],
                 Z     =inputs[3],
                 D     =inputs[4],
                 theta =inputs[5],
                 b1    =inputs[6],
                 b2    =inputs[0],
                 b3    =inputs[0],
                 d1    =21,
                 d2    =num_days+1,
                 d3    =num_days+1,
                 theta1=inputs[7],
                 theta2=inputs[5],
                 theta3=inputs[5])

       k = 0
       for c in cantons_:
          IU0[c] =   inputs[8+k]
          E0 [c] = 3*inputs[8+k]
          k += 1

       S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
       y0 = S0 + E0 + IR0 + IU0 + N0

       # Run the ODE solver.
       solver = libepidemics.cantons.seiin_interventions.Solver(data.to_cpp())
       y0 = libepidemics.cantons.seiin_interventions.State(y0)
       return solver.solve(params, y0, t_eval=range(1, num_days + 1))


    if L == 12 + ic_cantons: #case 3
       '''
       print("Running case 3 with:")
       print("beta  =",inputs[ 0])
       print("mu    =",inputs[ 1])
       print("alpha =",inputs[ 2])
       print("Z     =",inputs[ 3])
       print("D     =",inputs[ 4])
       print("theta =",inputs[ 5])
       print("b1    =",inputs[ 6])
       print("b2    =",inputs[ 7])
       print("d1    =",inputs[ 8])
       print("d2    =",inputs[ 9])
       print("theta1=",inputs[10])
       print("theta2=",inputs[11])
       '''

       params = libepidemics.cantons.seiin_interventions.Parameters(
                 beta  =inputs[0],
                 mu    =inputs[1],
                 alpha =inputs[2],
                 Z     =inputs[3],
                 D     =inputs[4],
                 theta =inputs[5],
                 b1    =inputs[6],
                 b2    =inputs[7],
                 b3    =inputs[0],
                 d1    =inputs[8],
                 d2    =inputs[9],
                 d3    =num_days+1,
                 theta1=inputs[10],
                 theta2=inputs[11],
                 theta3=inputs[5])

       k = 0
       for c in cantons_:
          IU0[c] =   inputs[12+k]
          E0 [c] = 3*inputs[12+k]
          k += 1

       S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
       y0 = S0 + E0 + IR0 + IU0 + N0

       # Run the ODE solver.
       solver = libepidemics.cantons.seiin_interventions.Solver(data.to_cpp())
       y0 = libepidemics.cantons.seiin_interventions.State(y0)
       return solver.solve(params, y0, t_eval=range(1, num_days + 1))


    params = libepidemics.cantons.seiin_interventions.Parameters(
              beta  =inputs[0],
              mu    =inputs[1],
              alpha =inputs[2],
              Z     =inputs[3],
              D     =inputs[4],
              theta =inputs[5],
              b1    =inputs[6],
              b2    =inputs[7],
              b3    =inputs[0],
              d1    =inputs[8],
              d2    =inputs[9],
              d3    =num_days+1,
              theta1=inputs[10],
              theta2=inputs[11],
              theta3=inputs[5])

    k = 0
    for c in cantons_:
       IU0[c] =   inputs[12+k]
       E0 [c] = 3*inputs[12+k]
       k += 1

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libepidemics.cantons.seiin_interventions.Solver(data.to_cpp())
    y0 = libepidemics.cantons.seiin_interventions.State(y0)
    return solver.solve(params, y0, t_eval=range(1, num_days + 1))
