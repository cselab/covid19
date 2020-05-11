import libepidemics

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx

def solve_sir(beta, gamma, y0, t_eval, *, N):
    """Solve the SIR equation with derivatives manually.

    Arguments:
        beta, gamma: model parameters
        y0: (S0, I0, R0) initial condition at t=0
        t_eval: A list of times `t` to return the values of the ODE at.
        N: population, model data

    Returns:
        A list of 9-tuples (S, I, R, dS/dbeta, dI/dbeta, dR/dbeta, dS/dgamma, dI/dgamma, dR/dgamma),
        one tuple for each element of t_eval.
    """
    def rhs(t, y):
        S, I, R, dSdbeta, dIdbeta, dRdbeta, dSdgamma, dIdgamma, dRdgamma = y

        A = beta * S * I / N
        B = gamma * I
        dAdbeta  = (S * I + beta * dSdbeta * I + beta * S * dIdbeta) / N
        dAdgamma = beta / N * (dSdgamma * I + S * dIdgamma)
        dBdbeta  = gamma * dIdbeta
        dBdgamma = I + gamma * dIdgamma

        dS = -A
        dI = A - B
        dR = B
        ddSdbeta = -dAdbeta
        ddIdbeta = dAdbeta - dBdbeta
        ddRdbeta = dBdbeta
        ddSdgamma = -dAdgamma
        ddIdgamma = dAdgamma - dBdgamma
        ddRdgamma = dBdgamma
        return (dS, dI, dR, ddSdbeta, ddIdbeta, ddRdbeta, ddSdgamma, ddIdgamma, ddRdgamma)

    # Initial derivatives are 0.
    results = solve_ivp(rhs, y0=(*y0, 0, 0, 0, 0, 0, 0),
                        t_span=(0, t_eval[-1]), t_eval=t_eval)
    results = np.transpose(results.y)
    assert len(results) == len(t_eval)
    assert len(results[0]) == 9
    return results



class TestCountrySIR(TestCaseEx):
    def test_sir(self):
        """Test the C++ autodiff implementation of the SIR model."""
        sir = libepidemics.country.sir
        data   = libepidemics.country.ModelData(N=1000000)
        solver = sir.Solver(data)
        params = sir.Parameters(beta=0.2, gamma=0.1)

        y0 = (10., 1., 2.)  # S, I, R.
        t_eval = [0, 0.3, 0.6, 1.0]
        initial = sir.State(y0)
        py_result = solve_sir(params.beta, params.gamma, y0=y0, t_eval=t_eval, N=data.N)
        cpp_result = solver.solve(params, initial, t_eval=t_eval)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, cpp in zip(py_result[1:], cpp_result[1:]):
            # See common.TestCaseEx.assertRelative
            self.assertRelative(py[0], cpp.S().val(), tolerance=1e-5)
            self.assertRelative(py[1], cpp.I().val(), tolerance=1e-5)
            self.assertRelative(py[2], cpp.R().val(), tolerance=1e-5)
            self.assertRelative(py[3], cpp.S().d(0),  tolerance=1e-5)
            self.assertRelative(py[4], cpp.I().d(0),  tolerance=1e-5)
            self.assertRelative(py[5], cpp.R().d(0),  tolerance=1e-5)
            self.assertRelative(py[6], cpp.S().d(1),  tolerance=1e-5)
            self.assertRelative(py[7], cpp.I().d(1),  tolerance=1e-5)
            self.assertRelative(py[8], cpp.R().d(1),  tolerance=1e-5)
