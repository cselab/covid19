from unittest import TestCase

import libepidemics
import numpy as np

def flatten(M: np.ndarray):
    """Flatten a numpy matrix."""
    return sum(M.tolist(), [])


def intervention_beta(t: float, p):
    """Compute the beta at time t, considering the intervetions from the
    parameter set p. Analoguous to the C++ function in intervention.h"""
    # Sanity check:
    p.tact
    p.dtact
    p.beta
    p.kbeta

    t0 = p.tact - 0.5 * p.dtact
    if t < t0:
        return p.beta
    elif t < p.tact + 0.5 * p.dtact:
        return (1. - (t - t0) / p.dtact * (1. - p.kbeta)) * p.beta
    else:
        return p.kbeta * p.beta


def gen_canton_design_parameters(K, days):
    """Return a random canton DesignParameters object.

    Arguments:
        K: number of cantons
        days: number of days (for commute data)
    """
    np.random.seed(12345)
    Mij = 1000 * np.random.rand(K, K)
    for i in range(K):
        Mij[i, i] = 0.0  # No self-flow.

    Cij = 1000 * np.random.rand(K, K)
    for i in range(K):
        Cij[i, i] = 0.0;

    return libepidemics.cantons.DesignParameters(
            region_keys=["C" + str(k) for k in range(K)],
            Ni=(1e6 + 1e6 * np.random.rand(K)).tolist(),
            Mij=flatten(Mij),
            Cij=flatten(Cij),
            ext_com_iu=flatten(100 * np.random.rand(days, K)),
            Ui=100 * np.random.rand(days))


class TestCaseEx(TestCase):
    def assertRelative(self, a, b, tolerance):
        if a == b:  # In the case both are 0.
            return
        relative = abs((a - b) / abs(a))
        if relative > tolerance:
            self.fail(f"assertRelative failed: |{a} - {b}| relative "
                      f"error {relative} larger than tolerance {tolerance}.")
