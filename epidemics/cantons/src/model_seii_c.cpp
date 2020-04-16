#include "model_seii_c.h"
#include "common.hh"

namespace seii_c {

void Solver::rhs([[maybe_unused]] int day, Parameters p, const State &x, State &dxdt) const
{
    const double * __restrict__ invNi = modelData_.invNi.data();
    for (size_t i = 0; i < modelData_.numRegions; ++i) {
        // double external_cases = modelData_.getExternalCasesAt(day, i);
        double sumIC_N = 0.0;
        for (size_t j = 0; j < modelData_.numRegions; ++j)
            sumIC_N += x.Iu(j) * this->C(i, j) * invNi[j];
        double sum_SC_N = 0.0;
        for (size_t j = 0; j < modelData_.numRegions; ++j)
            sum_SC_N += x.S(j) * this->C(j, i) * invNi[j];
        // printf("i=%zu invNi=%lg sumIC_N=%lg sum_SC_N=%lg\n", i, invNi[i], sumIC_N, sum_SC_N);
        double A = p.beta * invNi[i] * (
                x.S(i) * (x.Iu(i) + p.nu * sumIC_N)
                + p.nu * x.Iu(i) * sum_SC_N);
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
