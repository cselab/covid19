import libepidemics

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx

def solve_sir(beta, gamma, tact, dtact, kbeta, y0, t_eval, *, N):
    """Solve the SIR equation with interventions.

    Arguments:

    Returns:
    
    """
    def rhs(t, y):
        S, I, R = y

        betareal = None
        if (t < tact):
            betareal = beta
        elif (t < tact + dtact):
            betareal = (1. - (t - tact) / dtact * (1. - kbeta)) * beta
        else:
            betareal = kbeta*beta


        A = betareal * S * I / N
        B = gamma * I

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
        sir_int = libepidemics.country.sir_int
        data    = libepidemics.country.ModelData(N=100500)
        solver  = sir_int.Solver(data)
        params  = sir_int.Parameters(beta=0.9, gamma=0.6, tact=4.0, dtact=2.0, kbeta=0.5)

        y0 = (1e5, 1., 200.)  # S, I, R.
        t_eval    = [0, 0.3, 0.6, 1.0, 5.0, 10.0, 20.0]
        initial   = sir_int.State(y0)
        py_result = solve_sir(params.beta, params.gamma, params.tact, params.dtact, params.kbeta, y0=y0, t_eval=t_eval, N=data.N)
        cpp_result_noad = solver.solve   (params, initial, t_eval=t_eval, dt=0.01)
        cpp_result_ad   = solver.solve_ad(params, initial, t_eval=t_eval, dt=0.01)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, noad, ad in zip(py_result[1:], cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            self.assertRelative(noad.S(), ad.S().val(), tolerance=1e-12)
            self.assertRelative(noad.I(), ad.I().val(), tolerance=1e-12)
            self.assertRelative(noad.R(), ad.R().val(), tolerance=1e-12)
            self.assertRelative(py[0], ad.S().val(), tolerance=1e-6)
            self.assertRelative(py[1], ad.I().val(), tolerance=1e-6)
            self.assertRelative(py[2], ad.R().val(), tolerance=1e-6)
