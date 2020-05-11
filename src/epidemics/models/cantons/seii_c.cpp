#include "seii_c.h"
#include "common.hh"

namespace epidemics {
namespace cantons {
namespace seii_c {

void Solver::rhs([[maybe_unused]] int day, Parameters p, const State &x, State &dxdt) const
{
    const double * __restrict__ invNi = modelData_.invNi.data();
    for (size_t i = 0; i < modelData_.numRegions; ++i) {
        double sumIC_N = modelData_.getExternalCommutersIu(day, i);
        for (size_t j = 0; j < modelData_.numRegions; ++j)
            sumIC_N += x.Iu(j) * this->C_plus_Ct(i, j) * invNi[j];
        // printf("i=%zu invNi=%lg sumIC_N=%lg sum_SC_N=%lg\n", i, invNi[i], sumIC_N, sum_SC_N);
        double A = p.beta * x.S(i) * invNi[i] * (x.Iu(i) + p.nu * sumIC_N);
        double E_Z = x.E(i) / p.Z;

        double dS = -A;
        double dE = A - E_Z;
        double dIr = p.alpha * E_Z - x.Ir(i) / p.D;
        double dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;

        dxdt.S(i) = dS;
        dxdt.E(i) = dE;
        dxdt.Ir(i) = dIr;
        dxdt.Iu(i) = dIu;
    }
}

}  // namespace seii_c
}  // namespace cantons
}  // namespace epidemics
