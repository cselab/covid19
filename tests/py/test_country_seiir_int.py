import libepidemics

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx, intervention_beta

def solve_seiir(p, y0, t_eval, *, N):
    """Solve the SEIIR equation with interventions.

    Arguments:

    Returns:
    
    """
    def rhs(t, y):
        S, E, Ir, Iu, R = y
        
        betareal = intervention_beta(t, p)

        C1 = betareal * S * Ir / N
        C2 = betareal * S * (p.mu * Iu) / N
        C3 = p.alpha * E / p.Z
        C4 = (1 - p.alpha) * E / p.Z

        dS  = -C1 - C2
        dE  = C1 + C2 - E / p.Z
        dIr = C3 - Ir / p.D
        dIu = C4 - Iu / p.D
        dR = -dS - dE - dIr - dIu

        return (dS, dE, dIr, dIu, dR)


    results = solve_ivp(rhs, y0=y0, t_span=(0, t_eval[-1]), t_eval=t_eval, rtol=1e-9, atol=1e-9)
    results = np.transpose(results.y)
    assert len(results) == len(t_eval)
    assert len(results[0]) == 5
    return results



class TestCountrySIR(TestCaseEx):
    def test_sir(self):
        """Test the C++ autodiff implementation of the SIR model with interventions."""
        seiir_int  = libepidemics.country.seiir_int
        data       = libepidemics.country.ModelData(N=100500)
        solver     = seiir_int.Solver(data)
        params     = seiir_int.Parameters(beta=0.2, mu=0.1, alpha=0.15, Z=5.3, D=3.2, tact=10.0, dtact=2.0, kbeta=0.5)

        y0 = (1e5, 0.0, 1., 0.0, 200.)  # S, E, Ir, Iu, R
        t_eval    = [0, 0.3, 0.6, 1.0, 5.0, 10.0, 20.0]
        initial   = seiir_int.State(y0)
        py_result = solve_seiir(params, y0=y0, t_eval=t_eval, N=data.N)
        cpp_result_noad = solver.solve          (params, initial, t_eval=t_eval, dt=0.01)
        cpp_result_ad   = solver.solve_params_ad(params, initial, t_eval=t_eval, dt=0.01)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, noad, ad in zip(py_result[1:], cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            self.assertRelative(noad.S(), ad.S().val(), tolerance=1e-12)
            self.assertRelative(noad.E(), ad.E().val(), tolerance=1e-12)
            self.assertRelative(noad.Ir(), ad.Ir().val(), tolerance=1e-12)
            self.assertRelative(noad.Iu(), ad.Iu().val(), tolerance=1e-12)
            self.assertRelative(noad.R(), ad.R().val(), tolerance=1e-12)
            self.assertRelative(py[0], ad.S().val(), tolerance=1e-6)
            self.assertRelative(py[1], ad.E().val(), tolerance=1e-6)
            self.assertRelative(py[2], ad.Ir().val(), tolerance=1e-6)
            self.assertRelative(py[3], ad.Iu().val(), tolerance=1e-6)
            self.assertRelative(py[4], ad.R().val(), tolerance=1e-6)
