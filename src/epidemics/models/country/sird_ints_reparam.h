#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace sird_ints_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 5;

    T R0;     /// Transmission rate.
    T D;      /// Recovery time.
    T eps;    /// Mortality rate.
    T tact;   /// Day of intervention.
    T kbeta;  /// Multiplicator beta after intervention.
};

/// SIRD state has 4 elements: S, I, R.
template <typename T>
struct State : StateBase<T, 4> {
    using StateBase<T, 4>::StateBase;

    T &S() { return this->v_[0]; }
    T &I() { return this->v_[1]; }
    T &R() { return this->v_[2]; }
    T &D() { return this->v_[3]; }

    const T &S() const { return this->v_[0]; }
    const T &I() const { return this->v_[1]; }
    const T &R() const { return this->v_[2]; }
    const T &D() const { return this->v_[3]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {
        double invN = 1. / dp_.N;
        auto invD = 1. / p.D;

        T r0 = intervention_step(p.R0, t, p.kbeta, p.tact);
        
        auto A = invN * r0 * invD * x.I() * x.S();
        auto B = invD * x.I();

        dxdt.S() = -A;
        dxdt.I() = A - B;
        dxdt.R() = (1.0-p.eps)*B;
        dxdt.D() = p.eps*B;
    }
};

}  // namespace sird_ints_reparam
}  // namespace country
}  // namespace epidemics
