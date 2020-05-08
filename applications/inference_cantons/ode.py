from scipy.integrate import solve_ivp
import numpy as np

# also appends `sys.path` by `build/`
from epidemics.cantons.py.model import ModelData
import libsolver


def smooth_trans(u0, u1, t, tc, teps):
    """
    Smooth transition from u0 to u1 in interval `tc - teps < t < tc + teps`.
    """
    t0 = tc - teps
    t1 = tc + teps
    return u0 if t <= t0 else u1 if t > t1 else \
        u0 + (u1 - u0) * (1 - np.cos(np.pi/(t1 - t0)*(t - t0))) * 0.5


class Ode:
    params_fixed = dict()
    params_prior = dict()

    def solve(self, params, t_span, y0, t_eval):
        """
        Solves model equations.
        params: dict()
            Parameters.
        t_span: `2-tuple of floats`
            Interval of integration.
        y0: `array_like`, (n_vars, n_regions)
            Initial state.
        t_eval: `array_like`, (nt)
            Times at which to store the computed solution.
        Returns:
        y: `array_like`, (n_vars, n_regions, nt)
            Solution.
        """
        raise NotImplementedError()

    def solve_S_I_Icum(self, params, t_span, y0, t_eval):
        """
        Solves model equations converting the initial state from [S, I]
        and the solution to [S, I, Icum]
        S (susceptible), I (current infected), Icum (cumulative infected)

        params: dict()
            Parameters.
        t_span: `2-tuple of floats`
            Interval of integration.
        si0: `array_like`, (2, n_regions)
            Initial state.
        t_eval: `array_like`, (nt)
            Times at which to store the computed solution.
        Returns:
        s_i_icum: `array_like`, (3, n_regions, nt)
            Solution.
        """
        raise NotImplementedError()


class Sir(Ode):
    params_fixed = {'R0': 1.002, 'gamma': 60.}
    params_prior = {'R0': (0.5, 3), 'gamma': (1, 100)}

    def solve(self, params, t_span, y0, t_eval):
        """
        y0: `array_like`, (2, n_regions)
        """
        def rhs(t, y_flat):
            S, I = y_flat.reshape(2, -1)
            N = params['N']
            gamma = params['gamma']
            beta = params['R0'] * gamma
            c1 = beta * S * I / N
            c2 = gamma * I
            dSdt = -c1
            dIdt = c1 - c2
            return np.array((dSdt, dIdt)).flatten()

        y0 = np.array(y0)
        n_vars = 2
        assert y0.shape[0] == n_vars
        n_regions = y0.shape[1]
        y0_flat = y0.flatten()
        sol = solve_ivp(rhs, t_span=t_span, y0=y0_flat, t_eval=t_eval)
        y_flat = sol.y
        y = y_flat.reshape(n_vars, n_regions, len(t_eval))
        return y

    def solve_S_I_Icum(self, params, t_span, si0, t_eval):
        S, I = self.solve(params, t_span, si0, t_eval)
        N = params['N']
        Icum = N[:, None] - S
        return np.array((S, I, Icum))


class Seir(Ode):
    params_fixed = {
        'R0': 1.33,
        'Z': 1.5,
        'D': 2,
        'nu': 0.65,
        'theta_a': 0,
        'theta_b': 0.06,
        'tact': 31.,
        'kbeta': 0.5,
    }
    params_prior = {
        'R0': (0.5, 4),
        'Z': (0.1, 10),
        'D': (0.1, 10),
        'nu': (0.1, 10.),
        'theta_a': (0., 0.1),
        'theta_b': (0., 0.1),
        'tact': (0, 60),
        'kbeta': (0., 1.),
    }
    def __init__(self):
        for i in range(26):
            var = "beta_corr{:}".format(i)
            self.params_fixed[var] = 0
            self.params_prior[var] = (-0.5, 0.5)

    def solve(self, params, t_span, y0, t_eval):
        def rhs(t, y_flat):
            S, E, I = y_flat.reshape(3, -1)
            N = params['N']
            C = params['C']
            Z = params['Z']
            D = params['D']
            beta = params['R0'] / D
            beta = smooth_trans(beta, beta * params['kbeta'], t,
                                params['tact'], 0)
            nu = params['nu']
            k1 = I + nu * np.dot(C + C.T, I / N)
            A = beta * S * k1 / N
            #A = beta * k1 # XXX linearized
            E_Z = E / Z
            dS = -A
            dE = A - E_Z
            dI = E_Z - I / D
            return np.array((dS, dE, dI)).flatten()

        y0 = np.array(y0)
        n_vars = 3
        assert y0.shape[0] == n_vars
        n_regions = y0.shape[1]
        y0_flat = y0.flatten()
        sol = solve_ivp(rhs, t_span=t_span, y0=y0_flat, t_eval=t_eval)
        y_flat = sol.y
        y = y_flat.reshape(n_vars, n_regions, len(t_eval))
        return y

    def solve_S_I_Icum(self, params, t_span, si0, t_eval):
        S0, I0 = np.array(si0)
        E0 = np.zeros_like(S0)
        S, E, I = self.solve(params, t_span, [S0, E0, I0], t_eval)
        N = params['N']
        Icum = N[:, None] - S - E
        return np.array((S, I, Icum))


class SeirCpp(Seir):
    def solve(self, params, t_span, y0, t_eval):
        beta = params['R0'] / params['D']
        N = params['N']
        keys = list(map(str, range(len(N))))
        Cij = np.array(params['C'])
        Mij = np.zeros_like(Cij)

        y0 = np.array(y0).astype(float)
        n_vars = 3
        assert y0.shape[0] == n_vars
        n_regions = y0.shape[1]
        n_days = int(max(t_span)) + 1

        data = ModelData(keys, N, Mij, Cij)
        src = np.zeros(n_regions)
        src += params['theta_a'] * np.array(params['Qa'])
        src += params['theta_b'] * np.array(params['Qb'])
        data.ext_com_Iu = [src] * n_days
        data.Ui = [0] * n_regions
        for var in ['beta_corr' + str(i) for i in range(26)]:
            for i in params["beta_corr_regions"].get(var, []):
                data.Ui[i] = params[var]

        solver = libsolver.solvers.sei_c.Solver(data.to_cpp())

        p = libsolver.solvers.sei_c.Parameters(beta=beta,
                                               nu=params['nu'],
                                               Z=params['Z'],
                                               D=params['D'],
                                               tact=params['tact'],
                                               kbeta=params['kbeta'])
        sol = solver.solve(p, list(y0.flatten()), n_days)

        y = np.zeros((n_vars, n_regions, len(t_eval)))

        for i, t in enumerate(t_eval):
            q = min(len(sol) - 1, int(len(sol) * t / n_days))
            y[0, :, i] = sol[q].S()
            y[1, :, i] = sol[q].E()
            y[2, :, i] = sol[q].I()

        return y
