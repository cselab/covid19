import libepidemics

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx

def solve_seir(params, y0, t_eval, *, N):
    """Solve the SIR equation with derivatives manually.

    Arguments:
        params: a tuple (beta, gamma, a )
        y0: (S0, E0, I0, R0) initial condition at t=0
        t_eval: A list of times `t` to return the values of the ODE at.
        N: population, model data

    Returns:
        A list of 4-tuples (S, E, I, R),
        one tuple for each element of t_eval.
    """
    beta, gamma, a = params
    def rhs(t, y):
        S, E, I, R = y

        A = beta * S * I / N
        B = a * E
        C = gamma * I

        dS = -A
        dE = A - B
        dI = B - C
        dR = C
        return (dS, dE, dI, dR)

    results = solve_ivp(rhs, y0=y0, t_span=(0, t_eval[-1]), t_eval=t_eval, rtol=1e-9, atol=1e-9)
    results = np.transpose(results.y)
    assert len(results) == len(t_eval)
    assert len(results[0]) == 4
    return results



class TestCountrySIR(TestCaseEx):
    def test_sir(self):
        """Test the C++ autodiff implementation of the SIR model."""
        seir   = libepidemics.country.seir
        data   = libepidemics.country.ModelData(N=100500)
        solver = seir.Solver(data)
        params = seir.Parameters(beta=0.2, gamma=0.1, a = 1.0/3.0)

        y0 = (1e5, 0.0, 1., 200.)  # S, E, I, R.
        t_eval = [0, 0.3, 0.6, 1.0, 5.0, 10.0, 20.0]
        initial = seir.State(y0)
        py_result = solve_seir(params, y0=y0, t_eval=t_eval, N=data.N)
        cpp_result_noad = solver.solve          (params, initial, t_eval=t_eval, dt=0.1)
        cpp_result_ad   = solver.solve_params_ad(params, initial, t_eval=t_eval, dt=0.1)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, noad, ad in zip(py_result[1:], cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            self.assertRelative(noad.S(), ad.S().val(), tolerance=1e-12)
            self.assertRelative(noad.E(), ad.E().val(), tolerance=1e-12)
            self.assertRelative(noad.I(), ad.I().val(), tolerance=1e-12)
            self.assertRelative(noad.R(), ad.R().val(), tolerance=1e-12)
            self.assertRelative(py[0], ad.S().val(), tolerance=1e-7)
            self.assertRelative(py[1], ad.E().val(), tolerance=1e-7)
            self.assertRelative(py[2], ad.I().val(), tolerance=1e-7)
            self.assertRelative(py[3], ad.R().val(), tolerance=1e-7)
