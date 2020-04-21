#include "model_seiin.h"
#include "common.hh"

namespace seiin {

void Solver::rhs(int day, Parameters p, const State &x, State &dxdt) const
{
    for (size_t i = 0; i < modelData_.numRegions; ++i) {
        double extComIu = modelData_.getExternalCommutersIu(day, i);
        double A = p.beta * x.S(i) / x.N(i) * (x.Ir(i) + extComIu);
        double B = p.beta * x.S(i) / x.N(i) * p.mu * x.Iu(i);
        double E_Z = x.E(i) / p.Z;

        double dS = -(A + B);
        double dE = A + B - E_Z;
        double dIr = p.alpha * E_Z - x.Ir(i) / p.D;
        double dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;
        double dN = 0;

        double inv = 1 / (x.N(i) - x.Ir(i));
        //for (size_t j = 0; j < modelData_.numRegions; ++j) { // XXX
        for (size_t j : this->nonzero_Mij(i)) {
            double Tij = this->M(i, j) / (x.N(j) - x.Ir(j));
            double Tji = this->M(j, i) * inv;
            dS += p.theta * (Tij * x.S(j) - Tji * x.S(i));
            dE += p.theta * (Tij * x.E(j) - Tji * x.E(i));
            // Documented infected people are in quarantine, they do not move around.
            // dIr += p.theta * (Tij * x.Ir(j) - Tji * x.Ir(i));
            dIu += p.theta * (Tij * x.Iu(j) - Tji * x.Iu(i));
            dN += p.theta * (this->M(i, j) - this->M(j, i));
        }

        dxdt.S(i) = dS;
        dxdt.E(i) = dE;
        dxdt.Ir(i) = dIr;
        dxdt.Iu(i) = dIu;
        dxdt.N(i) = dN;
    }
}

}  // namespace seiin
