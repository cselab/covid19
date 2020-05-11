#include "base.hh"
#include "sir.h"

namespace epidemics {
namespace country {
namespace sir {

void Solver::rhs(double /* t */,
                 Parameters p,
                 const State &x,
                 State &dxdt) const
{
    auto beta  = make_ad(p.beta,  1, 0);
    auto gamma = make_ad(p.gamma, 0, 1);

    double invN = 1. / data_.N;
    auto A = invN * beta * x.I() * x.S();
    auto B = gamma * x.I();

    dxdt.S() = -A;
    dxdt.I() = A - B;
    dxdt.R() = B;
}

}  // namespace sir
}  // namespace country
}  // namespace epidemics
