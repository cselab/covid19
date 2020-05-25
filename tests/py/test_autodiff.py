from common import TestCaseEx, gen_canton_model_data

import libepidemics
from libepidemics.country import sir
from libepidemics.cantons import seiin

from epidemics.tools.autodiff import country_custom_derivatives, cantons_custom_derivatives

class TestAutoDiff(TestCaseEx):
    def test_country_custom_derivatives(self):
        from test_country_sir import solve_sir

        # We test derivatives wrt [nothing, beta, gamma, S0, I0, R0, beta+gamma+S0+I0+R0].
        # One column == one derivative.
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
        cpp_results, cpp_der_results = country_custom_derivatives(
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

    def test_cantons_custom_derivatives(self):
        K = 3  # Number of cantons.
        data = gen_canton_model_data(K=K, days=0)
        solver = seiin.Solver(data)
        params = seiin.Parameters(beta=0.3, mu=0.7, alpha=0.03, Z=4.0, D=5.0, theta=0.789)

        # S..., E..., Ir..., Iu..., N....
        y0 = seiin.State([1.0e6, 0.9e6, 0.8e6, 1, 2, 3, 5, 6, 7, 0, 1, 2, 3000000, 2000000, 1000000])
        t_eval = [0., 0.3, 0.6, 1.]

        # Test derivatives wrt (theta, D, Z, alpha, mu, beta, theta+D+Z+alpha+mu+beta).
        # One column == one derivative.
        params_der = (
            [0, 0, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 1, 0, 1],
            [0, 0, 0, 1, 0, 0, 1],
            [0, 0, 1, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0, 1],
        )
        # We wouldn't have anything to compare to, so we don't test derivatives
        # wrt initial conditions.
        y0_der = ([0, 0, 0, 0, 0, 0, 0], ) * len(y0)

        static_ad = solver.solve_params_ad(params, y0, t_eval=t_eval, dt=0.1)
        custom_results, custom_der_results = cantons_custom_derivatives(
                solver, params, y0, params_der, y0_der, t_eval=t_eval, dt=0.1)

        self.assertEqual(custom_results.shape,     (len(t_eval), 5, K))
        self.assertEqual(custom_der_results.shape, (len(t_eval), 5, K, 7))

        def compare(_static, _val, _der):
            """Compare value and derivatives for one canton, one day, one state variable."""
            try:
                self.assertEqual(_static.val(), _val)
                self.assertEqual(_static.d(0), _der[5])
                self.assertEqual(_static.d(1), _der[4])
                self.assertEqual(_static.d(2), _der[3])
                self.assertEqual(_static.d(3), _der[2])
                self.assertEqual(_static.d(4), _der[1])
                self.assertEqual(_static.d(5), _der[0])
                self.assertRelative(_der[6], sum(_der[0:6]), tolerance=1e-12)
            except:
                print(f"static={_static}")
                print(f"custom={_val}; {list(_der)}")
                raise

        for static, custom, custom_der in zip(static_ad, custom_results, custom_der_results):
            for k in range(K):
                compare(static.S(k),  custom[0, k], custom_der[0, k])
                compare(static.E(k),  custom[1, k], custom_der[1, k])
                compare(static.Ir(k), custom[2, k], custom_der[2, k])
                compare(static.Iu(k), custom[3, k], custom_der[3, k])
                compare(static.N(k),  custom[4, k], custom_der[4, k])
