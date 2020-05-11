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

    double inv_N = 1. / data_.N;
    auto A = inv_N * beta * x.I() * x.S();
    auto B = gamma * x.I();

    dxdt.S() = -A;
    dxdt.I() = A - B;
    dxdt.R() = B;
}

}  // namespace sir
}  // namespace country
}  // namespace epidemics
