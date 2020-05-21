from common import TestCaseEx

import libepidemics
from libepidemics.country import sir

from epidemics.tools.autodiff import compute_custom_derivatives

class TestAutoDiff(TestCaseEx):
    def test_compute_custom_derivatives(self):
        from test_country_sir import solve_sir

        # We test derivatives wrt [nothing, beta, gamma, S0, I0, R0, beta+gamma+S0+I0+R0].
        params = (0.3, 0.1)
        y0 = (1000000, 12, 3)
        N = sum(y0)
        params_derivatives = (
            [0, 1, 0, 0, 0, 0, 1],
            [0, 0, 1, 0, 0, 0, 1],
        )
        y0_derivatives = (
            [0, 0, 0, 1, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 1],
            [0, 0, 0, 0, 0, 1, 1],
        )
        t_eval = [float(x) for x in range(20)]

        data = libepidemics.country.ModelData(N=N)
        solver = sir.Solver(data)

        py_results = solve_sir(params, y0, t_eval=t_eval, N=N)
        cpp_results, cpp_der_results = compute_custom_derivatives(
                solver, params, y0, params_derivatives, y0_derivatives, t_eval=t_eval, dt=0.1)

        self.assertEqual(cpp_results.shape,     (len(t_eval), solver.state_size()))
        self.assertEqual(cpp_der_results.shape, (len(t_eval), solver.state_size(), 7))


        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for py, cpp, cpp_der in zip(py_results, cpp_results, cpp_der_results):
            self.assertRelative(py[0], cpp[0], tolerance=1e-7)  # S
            self.assertRelative(py[1], cpp[1], tolerance=1e-7)  # I
            self.assertRelative(py[2], cpp[2], tolerance=1e-7)  # R

            self.assertEqual(list(cpp_der[:, 0]), [0, 0, 0])  # Derivative wrt nothing.
            self.assertRelative(py[3], cpp_der[0, 1], tolerance=1e-7)  # dS/dbeta
            self.assertRelative(py[4], cpp_der[1, 1], tolerance=1e-7)  # dI/dbeta
            self.assertRelative(py[5], cpp_der[2, 1], tolerance=1e-7)  # dR/dbeta
            self.assertRelative(py[6], cpp_der[0, 2], tolerance=1e-7)  # dS/dgamma
            self.assertRelative(py[7], cpp_der[1, 2], tolerance=1e-7)  # dI/dgamma
            self.assertRelative(py[8], cpp_der[2, 2], tolerance=1e-7)  # dR/dgamma

            self.assertRelative(1, cpp_der[0, 3], tolerance=1e-2)  # dS/dS0
            self.assertLess(abs(0 - cpp_der[1, 3]), 1e-2)  # dI/dS0
            self.assertLess(abs(0 - cpp_der[2, 3]), 1e-2)  # dR/dS0

            # Testing dS/dI0, dI/dI0, dR/dI0 is not that simple...

            self.assertEqual(0, cpp_der[0, 5])  # dS/dR0
            self.assertEqual(0, cpp_der[1, 5])  # dI/dR0
            self.assertEqual(1, cpp_der[2, 5])  # dR/dR0

            self.assertRelative(cpp_der[0, 6], cpp_der[0, 0:6].sum(), 1e-12)  # dS/d(everthing)
            self.assertRelative(cpp_der[1, 6], cpp_der[1, 0:6].sum(), 1e-12)  # dI/d(everthing)
            self.assertRelative(cpp_der[2, 6], cpp_der[2, 0:6].sum(), 1e-12)  # dR/d(everthing)
