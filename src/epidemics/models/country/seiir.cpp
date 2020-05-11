#include "base.hh"
#include "seiir.h"

namespace epidemics {
namespace country {
namespace seiir {

void Solver::rhs(double /* t */,
                 Parameters p,
                 const State &x,
                 State & __restrict__ dxdt) const
{
    auto beta  = make_ad(p.beta,  1, 0, 0, 0, 0);
    auto mu    = make_ad(p.mu,    0, 1, 0, 0, 0);
    auto alpha = make_ad(p.alpha, 0, 0, 1, 0, 0);
    auto Z     = make_ad(p.Z,     0, 0, 0, 1, 0);
    auto D     = make_ad(p.D,     0, 0, 0, 0, 1);
    auto invZ  = Z.inv();  // 1 / Z.
    auto invD  = D.inv();  // 1 / D.

    double invN = 1. / data_.N;

    auto C1 = invN * beta * x.S() * x.Ir();
    auto C2 = invN * beta * x.S() * (mu * x.Iu());
    auto C3 = alpha * (invZ * x.E());
    auto C4 = (1 - alpha) * (invZ * x.E());

    dxdt.S()  = -(C1 + C2);
    dxdt.E()  = C1 + C2 - invZ * x.E();
    dxdt.Ir() = C3 - invD * x.Ir();
    dxdt.Iu() = C4 - invD * x.Iu();
    dxdt.R()  = -dxdt.S() - dxdt.E() - dxdt.Ir() - dxdt.Iu();
}

}  // namespace seiir
}  // namespace country
}  // namespace epidemics
