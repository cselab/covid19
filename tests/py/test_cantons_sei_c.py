from common import TestCaseEx, gen_canton_design_parameters

import libepidemics.cantons.sei_c as sei_c

# TODO: Compare with a Python implementation.

class TestCantonsSEI_C(TestCaseEx):
    def test_sei_c(self):
        """Test the C++ implementation of the SEI_C model."""
        K = 3  # Number of cantons.
        dp = gen_canton_design_parameters(K=K, days=10)
        solver = sei_c.Solver(dp)
        params = sei_c.Parameters(beta=0.3, nu=0.7, Z=0.03, D=4.0, tact=5.0, kbeta=0.789)  # Random.

        y0 = (10, 11, 12, 1, 2, 3, 5, 6, 7)  # (S..., E..., I...)
        t_eval = [0., 0.3, 0.6, 1., 1.5, 2., 3., 4.]
        y0 = sei_c.State(y0)
        cpp_result_noad = solver.solve          (params, y0, t_eval=t_eval)
        cpp_result_ad   = solver.solve_params_ad(params, y0, t_eval=t_eval)

        # Skip t=0 because relative error is undefined. Removing t=0 from t_eval does not work.
        for noad, ad in zip(cpp_result_noad[1:], cpp_result_ad[1:]):
            # See common.TestCaseEx.assertRelative
            for k in range(K):
                self.assertRelative(noad.S(k), ad.S(k).val(), tolerance=1e-12)
                self.assertRelative(noad.E(k), ad.E(k).val(), tolerance=1e-12)
                self.assertRelative(noad.I(k), ad.I(k).val(), tolerance=1e-12)
