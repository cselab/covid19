import numpy as np

def country_custom_derivatives(solver, params, y0, params_der, y0_der, t_eval, **kwargs):
    """Solve the ODE by customizing which derivatives are computed.

    This function can be used to compute derivatives with respect to initial conditions.

    Arguments:
        solver: a country model Solver object
        params: a tuple or a Parameters object
        y0: a tuple or a State object
        params_der: derivative matrix (*)
        y0_der: derivative matrix (*)
        t_eval: times at which to return the value of state variables
        kwargs: forwarded to solver.solve_ad

    Returns:
        A tuple of (results 2D matrix, derivatives 3D matrix). See the example for details.

    (*) The number of rows of params_der and y0_der must match the number of
    parameters and number of state variables, respectively. The columns
    determine different derivatives to be computed. Matrices params_dev and
    y0_dev must have an equal number of columns.

    Example:
        # Before. computing derivatives wrt parameters only.
        params = model.Parameters(beta=10, gamma=20, delta=30)
        y0 = (100000, 2, 1, 3)
        results = solver.solve_ad(params=params, y0=y0, t_eval=t_eval)
        for result in results:
            print(result.S().val())
            print(result.S().d(0))

        # After. Compute custom derivatives wrt 5 different variables: 3 parameters, y0[1] and y2[2].
        params = (10, 20, 30)
        # OR
        params = model.Parameters(beta=10, gamma=20, delta=30)
        y0 = (100000, 2, 1, 3)
        params_der = (
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
        )
        y0_der = (
            [0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0],
        )
        out, out_ad = country_custom_derivatives(solver, params, y0, params_der, y0_der, t_eval=t_eval)
        assert out.shape    == (len(t_eval), 4)
        assert out_ad.shape == (len(t_eval), 4, 5)
    """
    model = solver.model
    L = model.Parameters.NUM_PARAMETERS
    state_size = solver.state_size()
    if len(params) != L:
        raise TypeError(f"Expected {L} elements in `params`, got {len(params)}.")
    if len(params_der) != L:
        raise TypeError(f"Expected {L} elements in `params_der`, got {len(params_der)}.")
    if len(y0) != state_size:
        raise TypeError(f"Expected {state_size} elements in `y0`, got {len(y0)}.")
    if len(y0_der) != state_size:
        raise TypeError(f"Expected {state_size} elements in `y0_der`, got {len(y0_der)}.")

    K = len(params_der[0])  # Number of derivatives to compute.
    if not all(len(p) == K for p in params_der):
        raise TypeError("Subarray length in `params_der` is not consistent.")
    if not all(len(yi) == K for yi in y0_der):
        raise TypeError("Subarray length in `y0_der` is not consistent with `params_der`.")

    def extend_to_L(array):
        l = len(array)
        return array if l == L else list(array) + [0] * (L - l)

    results = np.full((len(t_eval), len(y0)), -1.0)
    results_der = np.full((len(t_eval), len(y0), K), -1.0)
    for k0 in range(0, K, L):
        k1 = min(k0 + L, K)
        tmp_params = model.ParametersAD(*(
                model.AD(param, extend_to_L(param_der[k0:k1]))
                for param, param_der in zip(params, params_der)))
        tmp_y0 = model.StateAD([
                model.AD(yi, extend_to_L(yi_der[k0:k1]))
                for yi, yi_der in zip(y0, y0_der)])

        partial = solver.solve_ad(params=tmp_params, y0=tmp_y0, t_eval=t_eval, **kwargs)
        for t, day in enumerate(partial):
            for i in range(len(y0)):
                results[t, i] = day(i).val()
                results_der[t, i, k0:k1] = day(i).d()[:k1 - k0]  # .d() returns a list of numbers.
    return results, results_der
