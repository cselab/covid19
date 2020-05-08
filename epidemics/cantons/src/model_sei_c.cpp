#include "model_sei_c.h"
#include "common.hh"

namespace sei_c {

void Solver::rhs(int day, Parameters p, const State &x, State &dxdt) const
{
    double beta = p.beta;
    // FIXME: rhs should take `double day` (same in lambda `rhs`)
    //        then replace `day + 1` by `day`
    if (day + 1 > p.tact) {
      beta *= p.kbeta;
    }
    const double * __restrict__ invNi = modelData_.invNi.data();
    for (size_t i = 0; i < modelData_.numRegions; ++i) {
        double sumIC_N = 0;
        for (size_t j = 0; j < modelData_.numRegions; ++j) {
            sumIC_N += x.I(j) * this->C_plus_Ct(i, j) * invNi[j];
        }
        const double ext = modelData_.getExternalCommutersIu(day, i);
        const double beta_i = beta * (1 + modelData_.Ui[i]);
        const double A = beta_i * x.S(i) * invNi[i] * (
                x.I(i) + p.nu * sumIC_N + ext);
        const double E_Z = x.E(i) / p.Z;

        const double dS = -A;
        const double dE = A - E_Z;
        const double dI = E_Z - x.I(i) / p.D;

        dxdt.S(i) = dS;
        dxdt.E(i) = dE;
        dxdt.I(i) = dI;
    }
}

}  // namespace sei_c
