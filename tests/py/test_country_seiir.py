import libepidemics

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx

def solve_seiir(p, y0, t_eval, *, N):
    """Solve the SEIIR equation.

    Arguments:
        p: model parameters
        y0: (S0, I0, R0) initial condition at t=0
        t_eval: A list of times `t` to return the values of the ODE at.
        N: population, model data

    Returns:
        A list of 5-tuples (S, E, Ir, Iu, R), one tuple for each element of t_eval.
    """
    def rhs(t, y):
        S, E, Ir, Iu, R = y

        C1 = p.beta * S * Ir / N
        C2 = p.beta * S * (p.mu * Iu) / N
        C3 = p.alpha * E / p.Z
        C4 = (1 - p.alpha) * E / p.Z

        dS  = -C1 - C2
        dE  = C1 + C2 - E / p.Z
        dIr = C3 - Ir / p.D
        dIu = C4 - Iu / p.D
        dR = -dS - dE - dIr - dIu

        return (dS, dE, dIr, dIu, dR)

    # Initial derivatives are 0.
    results = solve_ivp(rhs, y0=y0, t_span=(0, t_eval[-1]), t_eval=t_eval, rtol=1e-10, atol=1e-10)
    results = np.transpose(results.y)
    assert len(results) == len(t_eval)
    assert len(results[0]) == 5
    return results



class TestCountrySEIIR(TestCaseEx):
    def test_seiir(self):
        
        seiir  = libepidemics.country.seiir
        data   = libepidemics.country.ModelData(N=100050)
        solver = seiir.Solver(data)
        params = seiir.Parameters(beta=p[0], mu=p[1], alpha=p[2], Z=p[3], D=p[4])

        y0 = (1e5, 1., 2., 3., 5.)  # S, E, Ir, Iu, R.
        t_eval = [0, 0.3, 0.6, 1.0, 5.0, 10.0, 20.0]
        initial = seiir.State(y0)
        py_result = solve_seiir(params, y0=y0, t_eval=t_eval, N=data.N)
        cpp_result_noad = solver.solve   (params, initial, t_eval=t_eval, dt=0.1)
        cpp_result_ad   = solver.solve_ad(params, initial, t_eval=t_eval, dt=0.1)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for t, py, noad, ad in zip(t_eval[1:], py_result[1:], cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            try:
                self.assertRelative(noad.S(),  ad.S().val(),  tolerance=1e-12)
                self.assertRelative(noad.E(),  ad.E().val(),  tolerance=1e-12)
                self.assertRelative(noad.Ir(), ad.Ir().val(), tolerance=1e-12)
                self.assertRelative(noad.Iu(), ad.Iu().val(), tolerance=1e-12)
                self.assertRelative(noad.R(),  ad.R().val(),  tolerance=1e-12)
                self.assertRelative(py[0], ad.S().val(),  tolerance=1e-7)
                self.assertRelative(py[1], ad.E().val(),  tolerance=1e-7)
                self.assertRelative(py[2], ad.Ir().val(), tolerance=1e-7)
                self.assertRelative(py[3], ad.Iu().val(), tolerance=1e-7)
                self.assertRelative(py[4], ad.R().val(),  tolerance=1e-7)
            except AssertionError:
                print("Failed for t = ", t)
                print("    py   =", list(py))
                print("    noad =", noad.S(), noad.E(), noad.Ir(), noad.Iu(), noad.R())
                print("    ad   =", ad.S(), ad.E(), ad.Ir(), ad.Iu(), ad.R())
                raise
