#include "model_sei_c.h"
#include "common.hh"

namespace sei_c {

void Solver::rhs(int day, Parameters p, const State &x, State &dxdt) const
{
    const double * __restrict__ invNi = modelData_.invNi.data();
    for (size_t i = 0; i < modelData_.numRegions; ++i) {
        double sumIC_N = modelData_.getExternalCommutersIu(day, i);
        for (size_t j = 0; j < modelData_.numRegions; ++j) {
            sumIC_N += x.I(j) * this->C_plus_Ct(i, j) * invNi[j];
        }
        const double A = p.beta * x.S(i) * invNi[i] * (x.I(i) + p.nu * sumIC_N);
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
