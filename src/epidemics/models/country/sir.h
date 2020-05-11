#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace sir {

struct Parameters {
    double beta;   /// Transmission rate.
    double gamma;  /// Recovery rate.
};

/// Each of S, I and R are autodiff objects with 2 derivatives.
using Element = AutoDiff<double, 2>;

/// SIR state has 3 elements: S, I, R.
struct State : StateBase<Element, 3> {
    using StateBase<Element, 3>::StateBase;

    Element &S() { return v_[0]; }
    Element &I() { return v_[1]; }
    Element &R() { return v_[2]; }

    const Element &S() const { return v_[0]; }
    const Element &I() const { return v_[1]; }
    const Element &R() const { return v_[2]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    void rhs(double t, Parameters p, const State &x, State &dxdt) const;
};

}  // namespace sir
}  // namespace country
}  // namespace epidemics
