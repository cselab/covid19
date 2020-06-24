#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),'..',  '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__),'..',  '..', '..', 'build'))
import numpy as np
import matplotlib.pyplot as plt
from copy import copy

from epidemics.epidemics import EpidemicsBase
from epidemics.tools.tools import save_file
import libepidemics


class Par:
    R0 = 2.4
    gamma = 1. / 5.2
    tint = 30
    dint = 3
    kint = 0.2
    N = int(8e6)
    I0 = 100
    E0 = 100
    incubation = 5.
    tmax = 100


def solve_sir_r0(p):
    model = libepidemics.country.sir_int_r0
    data = libepidemics.country.ModelData(N=p.N)
    cppsolver = model.Solver(data)

    y0 = [p.N - p.I0, p.I0]

    params = model.Parameters(r0=p.R0,
                              gamma=p.gamma,
                              tact=p.tint,
                              dtact=p.dint,
                              kbeta=p.kint)

    s0, i0 = y0
    y0cpp = (s0, i0, 0.0)
    initial = model.State(y0cpp)
    t_eval = np.arange(p.tmax)
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt=0.1)
    S = np.zeros(len(cpp_res))
    for i, e in enumerate(cpp_res):
        S[i] = e.S()
    daily = -np.diff(S)
    return daily, copy(p)


def solve_sir(p):
    model = libepidemics.country.sir_int
    data = libepidemics.country.ModelData(N=p.N)
    cppsolver = model.Solver(data)

    y0 = [p.N - p.I0, p.I0]

    params = model.Parameters(beta=p.R0 * p.gamma,
                              gamma=p.gamma,
                              tact=p.tint,
                              dtact=p.dint,
                              kbeta=p.kint)
    S0 = p.N - p.I0
    y0cpp = (S0, p.I0, 0.0)
    initial = model.State(y0cpp)
    t_eval = np.arange(p.tmax)
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt=0.1)
    S = np.zeros(len(cpp_res))
    for i, e in enumerate(cpp_res):
        S[i] = e.S()
    daily = -np.diff(S)
    return daily, copy(p)


def solve_seir(p):
    model = libepidemics.country.seir_int
    data = libepidemics.country.ModelData(N=p.N)
    cppsolver = model.Solver(data)

    params = model.Parameters(beta=p.R0 * p.gamma,
                              a=1./p.incubation,
                              gamma=p.gamma,
                              tact=p.tint,
                              dtact=p.dint,
                              kbeta=p.kint)

    S0 = p.N - p.I0
    y0cpp = (S0, p.E0, p.I0, 0.0)
    initial = model.State(y0cpp)
    t_eval = np.arange(p.tmax)
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt=0.1)
    S = np.zeros(len(cpp_res))
    for i, e in enumerate(cpp_res):
        S[i] = e.S()
    daily = -np.diff(S)
    return daily, copy(p)

p = Par()
sir = solve_sir(p)
k = 2.3
p_seir = copy(p)
p_seir.R0 *= k
p_seir.kint /= k ** 2
seir = solve_seir(p_seir)

# R0seir from doc/linear_sir_vs_seir
R0seir = p.R0 * (1 + (p.R0 - 1) * p.gamma * p.incubation)
print("fitted R0 = ", p_seir.R0)
print("analytical R0 = ", R0seir)

#plt.plot(sir_r0, label="sirR0")
plt.plot(sir[0], label="sir,R0={:.3g},kin={:.3g}".format(sir[1].R0, sir[1].kint))
plt.plot(seir[0], label="seir,R0={:.3g},kin={:.3g}".format(seir[1].R0, seir[1].kint))
plt.yscale('log')
plt.legend()
plt.savefig("sir_vs_seir.pdf")

