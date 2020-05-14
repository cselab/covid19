import libepidemics.cantons.seiin as seiin

from scipy.integrate import solve_ivp
import numpy as np

from common import TestCaseEx, gen_canton_model_data

# TODO: Test external commuters! Set days to > 0.

def solve_seiin(params, y0, t_eval, *, data):
    """Solve the SEIIN equation.

    Arguments:
        params: model parameters
        y0: (S0, E0, Ir0, Iu0, N0) initial condition at t=0
        t_eval: A list of times `t` to return the values of the ODE at.
        Mij

    Returns:
        A list of 5-tuples, one tuple for each element of t_eval.
    """
    beta, mu, alpha, Z, D, theta = params
    K = data.num_regions
    Mij = np.array(data.Mij).reshape((K, K))

    colsumMij = np.sum(Mij, axis=0)
    rowsumMij = np.sum(Mij, axis=1)

    def rhs(t, y):
        assert len(y) == 5 * K, (len(y), K)
        S  = y[0 * K : 1 * K]
        E  = y[1 * K : 2 * K]
        Ir = y[2 * K : 3 * K]
        Iu = y[3 * K : 4 * K]
        N  = y[4 * K : 5 * K]

        # Here we try to avoid repeating any computation...
        tmpI = beta * S / N * (Ir + mu * Iu)
        tmpE_Z = E / Z
        tmpalphaE_Z = alpha * tmpE_Z

        dS = -tmpI
        dE = tmpI - tmpE_Z
        dIr = tmpalphaE_Z - Ir / D
        dIu = tmpE_Z - tmpalphaE_Z - Iu / D
        # dN = np.zeros(K)

        tmpNI = N - Ir
        tmpS_NI = S / tmpNI
        tmpE_NI = E / tmpNI
        tmpIu_NI = Iu / tmpNI
        dS  += theta * (np.matmul(Mij, tmpS_NI)  - colsumMij * tmpS_NI)
        dE  += theta * (np.matmul(Mij, tmpE_NI)  - colsumMij * tmpE_NI)
        dIu += theta * (np.matmul(Mij, tmpIu_NI) - colsumMij * tmpIu_NI)
        dN   = theta * (rowsumMij - colsumMij)  # This can be also set to 0.

        return [*dS, *dE, *dIr, *dIu, *dN]

    # Initial derivatives are 0.
    results = solve_ivp(rhs, y0=y0, t_span=(0, t_eval[-1]), t_eval=t_eval, rtol=1e-8, atol=1e-8)
    results = np.transpose(results.y)
    assert len(results) == len(t_eval)
    assert len(results[0]) == K * 5
    return results



class TestCantonsSEIIN(TestCaseEx):
    def test_seiin(self):
        """Test the C++ implementation of the SEIIN model."""
        K = 3  # Number of cantons.
        data = gen_canton_model_data(K=K, days=0)
        solver = seiin.Solver(data)
        params = seiin.Parameters(beta=0.3, mu=0.7, alpha=0.03, Z=4.0, D=5.0, theta=0.789)

        # S..., E..., Ir..., Iu..., N....
        y0 = (1.0e6, 0.9e6, 0.8e6, 1, 2, 3, 5, 6, 7, 0, 1, 2, 1000000, 2000000, 3000000)
        t_eval = [0., 0.3, 0.6, 1.]
        py_result = solve_seiin(params, y0=y0, t_eval=t_eval, data=data)
        y0 = seiin.State(y0)
        cpp_result_noad = solver.solve   (params, y0, t_eval=t_eval, dt=0.1)
        cpp_result_ad   = solver.solve_ad(params, y0, t_eval=t_eval, dt=0.1)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, noad, ad in zip(py_result[1:], cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            for k in range(K):
                self.assertRelative(noad.S(k),  ad.S(k).val(),  tolerance=1e-12)
                self.assertRelative(noad.E(k),  ad.E(k).val(),  tolerance=1e-12)
                self.assertRelative(noad.Ir(k), ad.Ir(k).val(), tolerance=1e-12)
                self.assertRelative(noad.Iu(k), ad.Iu(k).val(), tolerance=1e-12)
                self.assertRelative(noad.N(k),  ad.N(k).val(),  tolerance=1e-12)
                self.assertRelative(noad.S(k),  py[0 * K + k], tolerance=1e-7)
                self.assertRelative(noad.E(k),  py[1 * K + k], tolerance=1e-7)
                self.assertRelative(noad.Ir(k), py[2 * K + k], tolerance=1e-7)
                self.assertRelative(noad.Iu(k), py[3 * K + k], tolerance=1e-7)
                self.assertRelative(noad.N(k),  py[4 * K + k], tolerance=1e-7)
