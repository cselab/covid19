#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace seiir {

struct Parameters {
    double beta;   /// Transmission rate.
    double mu;     /// TODO
    double alpha;  /// TODO
    double Z;      /// TODO
    double D;      /// TODO
};

/// We want to differentiate each value wrt 5 parameters.
using Element = AutoDiff<double, 5>;

/// SEIIR state has 5 elements: S, E, I, I, R.
struct State : StateBase<Element, 5> {
    using StateBase<Element, 5>::StateBase;

    Element &S()  { return v_[0]; }
    Element &E()  { return v_[1]; }
    Element &Ir() { return v_[2]; }
    Element &Iu() { return v_[3]; }
    Element &R()  { return v_[4]; }

    const Element &S()  const { return v_[0]; }
    const Element &E()  const { return v_[1]; }
    const Element &Ir() const { return v_[2]; }
    const Element &Iu() const { return v_[3]; }
    const Element &R()  const { return v_[4]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    void rhs(double t, Parameters p, const State &x, State &dxdt) const;
};

}  // namespace seiir
}  // namespace country
}  // namespace epidemics

