#!/usr/bin/env python3

import torch

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.cantons.py.data import CANTON_POPULATION, fetch_canton_data
from epidemics.cantons.py.solver import Solver
from epidemics.cantons.py.plot_ode import plot_ode_results
import libsolver  # Must be AFTER .plot_ode (because of sys.path...).

CANTON_TO_INDEX, REFDATA = fetch_canton_data()

def get_params():
    """Get the 6 model parameters."""
    # Parameters: beta, mu, alpha, Z, D, theta
    # Li2020, Table 1
    return (1.12, 0., 1., 3.69, 3.47, 1.36)


def get_Mij():
    """Return the commute matrix."""
    Mij = 1000 * torch.rand((26, 26))
    Mij = Mij + Mij.t()  # No migrations.
    for i in range(26):
        Mij[i, i] = 0.0  # No self-flow.
    return Mij


def get_initial_state():
    """Prepare an initial state, a vector of 5*26 values."""
    # Common to any initial state.
    num_cantons = len(CANTON_TO_INDEX)
    S0 = [CANTON_POPULATION[canton] for canton in CANTON_TO_INDEX.keys()]
    E0 = [0] * num_cantons
    IR0 = [0] * num_cantons
    IU0 = [0] * num_cantons
    N0 = [CANTON_POPULATION[canton] for canton in CANTON_TO_INDEX.keys()]

    # 1 infected in Ticino.
    IR0[CANTON_TO_INDEX['TI']] = 1

    # The order must match the one in `Solver.solve.rhs`! Also, for the
    # visualization to work, the order must match the one in the C++ code!
    y0 = S0 + E0 + IR0 + IU0 + N0
    y0 = torch.FloatTensor(y0)
    return y0


def solve_and_visualize(y0, params, Mij, num_days):
    """Solve the ODE for the given set of parameters and initial values, and generate a movie."""
    solver = Solver(Mij)
    states = solver.solve(y0, params, num_days)

    # Convert to C++ State object, expected by `plot_ode_results`.
    states = [libsolver.State(state.tolist()) for state in states]
    plot_ode_results(states)


def main():
    y0 = get_initial_state()
    params = get_params()
    Mij = get_Mij()
    solve_and_visualize(y0, params, Mij, num_days=40)


if __name__ == '__main__':
    main()
