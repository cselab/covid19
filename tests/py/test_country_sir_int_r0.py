import libepidemics

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx

def solve_sir(p, y0, t_eval, *, N):
    """Solve the SIR equation with interventions.

    Arguments:

    Returns:
    
    """
    def rhs(t, y):
        S, I, R = y

        t0 = p.tact - 0.5 * p.dtact
        if t < t0:
            r0real = p.r0
        elif t < p.tact + 0.5 * p.dtact:
            r0real = (1. - (t - t0) / p.dtact * (1. - p.kbeta)) * p.r0
        else:
            r0real = p.kbeta * p.r0

        A = r0real * p.gamma * S * I / N
        B = p.gamma * I

        dS = -A
        dI = A - B
        dR = B
        return (dS, dI, dR)

    results = solve_ivp(rhs, y0=y0, t_span=(0, t_eval[-1]), t_eval=t_eval, rtol=1e-9, atol=1e-9)
    results = np.transpose(results.y)
    assert len(results) == len(t_eval)
    assert len(results[0]) == 3
    return results



class TestCountrySIR(TestCaseEx):
    def test_sir(self):
        """Test the C++ autodiff implementation of the SIR model with interventions."""
        sir_int_r0 = libepidemics.country.sir_int_r0
        dp         = libepidemics.country.DesignParameters(N=100500)
        solver     = sir_int_r0.Solver(dp)
        params     = sir_int_r0.Parameters(r0=0.99, gamma=0.6, tact=4.0, dtact=2.0, kbeta=0.5)

        y0 = (1e5, 1., 200.)  # S, I, R.
        t_eval    = [0, 0.3, 0.6, 1.0, 5.0, 10.0, 20.0]
        initial   = sir_int_r0.State(y0)
        py_result = solve_sir(params, y0=y0, t_eval=t_eval, N=dp.N)
        cpp_result_noad = solver.solve          (params, initial, t_eval=t_eval, dt=0.01)
        cpp_result_ad   = solver.solve_params_ad(params, initial, t_eval=t_eval, dt=0.01)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, noad, ad in zip(py_result[1:], cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            self.assertRelative(noad.S(), ad.S().val(), tolerance=1e-12)
            self.assertRelative(noad.I(), ad.I().val(), tolerance=1e-12)
            self.assertRelative(noad.R(), ad.R().val(), tolerance=1e-12)
            self.assertRelative(py[0], ad.S().val(), tolerance=1e-6)
            self.assertRelative(py[1], ad.I().val(), tolerance=1e-6)
            self.assertRelative(py[2], ad.R().val(), tolerance=1e-6)
