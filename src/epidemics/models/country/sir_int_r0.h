#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace sir_int_r0 {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 5;

    T r0;   /// Transmission rate.
    T gamma;  /// Recovery rate.
    T tact;   /// Day of intervention.
    T dtact;  /// Duration of intervention.
    T kbeta;  /// Multiplicator beta after intervention.
};

/// SIR state has 3 elements: S, I, R.
template <typename T>
struct State : StateBase<T, 3> {
    using StateBase<T, 3>::StateBase;

    T &S() { return this->v_[0]; }
    T &I() { return this->v_[1]; }
    T &R() { return this->v_[2]; }

    const T &S() const { return this->v_[0]; }
    const T &I() const { return this->v_[1]; }
    const T &R() const { return this->v_[2]; }
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
        
        T r0 = intervention(p.r0, t, p.kbeta, p.tact, p.dtact);
        
        auto A = invN * r0 * p.gamma * x.I() * x.S();
        auto B = p.gamma * x.I();

        dxdt.S() = -A;
        dxdt.I() = A - B;
        dxdt.R() = B;
    }
};

}  // namespace sir_int_r0
}  // namespace country
}  // namespace epidemics
