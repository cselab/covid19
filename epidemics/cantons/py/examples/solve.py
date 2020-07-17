#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os
import sys
import collections

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..',
                             '..'))

# also appends `sys.path` by `build/`
from epidemics.cantons.py.model import PyDesignParameters

import libepidemics

params = libepidemics.cantons.sei_c.Parameters(
        beta=1.12, nu=0., Z=3.69, D=3.47, tact=3)

n_regions = 2
population = [100] * n_regions

N0 = np.array(population)
E0 = np.zeros(n_regions)
I0 = np.zeros(n_regions) + 1

Mij = np.zeros((n_regions, n_regions))
Cij = np.zeros((n_regions, n_regions))

dp = PyDesignParameters(['ZH', 'TI'], [100, 200], Mij, Cij)

src = np.zeros(dp.num_regions)
dp.ext_com_Iu = [src]

solver = libepidemics.cantons.sei_c.Solver(dp.to_cpp())

S0 = N0 - E0 - I0
y0 = libepidemics.cantons.sei_c.State(np.vstack((S0, E0, I0)).flatten())
num_days = 10

r = solver.solve(params, y0, t_eval=range(1, num_days + 1))

print(list([s.S()] for s in r))
