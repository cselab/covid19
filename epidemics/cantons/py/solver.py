import torch

VARS_PER_REGION = 5

def RK4_step(f, y, t, dt):
    k1 = dt * f(y,t)
    k2 = dt * f(y + 0.5 * k1, t + 0.5 * dt)
    k3 = dt * f(y + 0.5 * k2, t + 0.5 * dt)
    k4 = dt * f(y + k3, t+dt)
    return y + k1/6 + k2/3 + k3/3 + k4/6, t + dt

def FE_step(f, y, t, dt):
    return y + dt * f(y,t), t + dt


def integrate(rhs_func, y0, num_days, iterations_per_day):
    """Integrate dy/dt = rhs_func from t=0 to t=num_days with y(0) = y0.

    Returns a list of state vectors, one for each day.
    """

    out = [y0.clone()]

    y = y0.clone()
    t = 0
    dt = 1 / iterations_per_day
    for day in range(num_days):
        for it in range(iterations_per_day):
            y += dt * rhs_func(y, t)
            t += dt
            # y, t = RK4_step(rhs_func,y,t,dt)
        out.append(y.clone())

    return out


class Solver:
    """PyTorch-based solver for the multi-region SEII model.

    Parameters of the model:
        Mij: A commute matrix.
        beta, mu, alpha, Z, D, theta: Various scalars.

    Implementation notice:
        There are 5 values integrated for each region/canton: S, E, Ir, Iu, N.
        We store the state as a concatenated vector [S..., E..., Ir..., Iu..., N...].

    Reference:
        "Substantial undocumented infection facilitates the rapid dissemination of novel coronavirus (SARS-CoV2)"
        Li et al., 20020, Supplementary Text, Page 13,
        https://science.sciencemag.org/content/early/2020/03/24/science.abb3221/
    """

    def __init__(self, Mij):
        assert len(Mij.size()) == 2
        assert Mij.size()[0] == Mij.size()[1]

        self.num_cantons = Mij.size()[0]
        self.Mij = Mij

    def solve(self, y0, params, num_days):
        assert len(y0) == VARS_PER_REGION * self.num_cantons

        num_cantons = self.num_cantons
        Mij = self.Mij
        colsumMij = torch.sum(Mij, dim=0)
        rowsumMij = torch.sum(Mij, dim=1)

        beta, mu, alpha, Z, D, theta = params

        def rhs(y, t):
            S  = y[0 * num_cantons : 1 * num_cantons]
            E  = y[1 * num_cantons : 2 * num_cantons]
            Ir = y[2 * num_cantons : 3 * num_cantons]
            Iu = y[3 * num_cantons : 4 * num_cantons]
            N  = y[4 * num_cantons : 5 * num_cantons]

            # Here we try to avoid repeating any computation...
            tmpI = beta * S / N * (Ir + mu * Iu)
            tmpE_Z = E / Z
            tmpalphaE_Z = alpha * tmpE_Z

            dS = -tmpI
            dE = tmpI - tmpE_Z
            dIr = tmpalphaE_Z - Ir / D
            dIu = tmpE_Z - tmpalphaE_Z - Iu / D
            # dN = torch.zeros(num_cantons)

            tmpNI = N - Ir
            tmpS_NI = S / tmpNI
            tmpE_NI = E / tmpNI
            tmpIu_NI = Iu / tmpNI
            dS  += theta * (torch.mv(Mij, tmpS_NI)  - colsumMij * tmpS_NI)
            dE  += theta * (torch.mv(Mij, tmpE_NI)  - colsumMij * tmpE_NI)
            dIu += theta * (torch.mv(Mij, tmpIu_NI) - colsumMij * tmpIu_NI)
            dN   = theta * (colsumMij - rowsumMij)  # This can be also set to 0.

            return torch.cat((dS, dE, dIr, dIu, dN))

        return integrate(rhs, y0, num_days, iterations_per_day=10)
