#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os
import sys
import collections

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))

# also appends `sys.path` by `build/`
from epidemics.cantons.py.model import ModelData

import libsolver

params = libsolver.solvers.seiin.Parameters(beta=1.12,
                                            mu=0.,
                                            alpha=1.,
                                            Z=3.69,
                                            D=3.47,
                                            theta=1.36)

n_regions = 2
population = [100] * n_regions

N0 = np.array(population)
E0 = np.zeros(n_regions)
IR0 = np.zeros(n_regions)
IU0 = np.zeros(n_regions)

Mij = np.zeros((n_regions, n_regions))
Cij = np.zeros((n_regions, n_regions))

data = ModelData(['ZH', 'TI'], [100, 200], Mij, Cij)

src = np.zeros(data.num_regions)
data.ext_com_Iu = [src]

solver = libsolver.solvers.seiin.Solver(data.to_cpp(), verbose=True)

S0 = N0 - E0 - IR0 - IU0
y0 = list(np.vstack((S0, E0, IR0, IU0, N0)).flatten())
num_days = 10

r = solver.solve(params, y0, num_days)

print(r[0].S())
