#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace sir_intexp_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 4;

    T R0;     /// Transmission rate.
    T D;      /// Recovery time.
    T tact;   /// Day of intervention.
    T k;      /// Decay.
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
        double invN = 1. / data_.N;
        auto invD = 1. / p.D;

        T r0;
        if (t <= p.tact) {
           r0 = p.R0;
        }
        else {
           using std::exp;
           double dt = p.tact - t; // careful: diff will not work
           r0 = p.R0*exp(p.k*dt);
        }
        auto A = invN * r0 * invD * x.I() * x.S();
        auto B = invD * x.I();

        dxdt.S() = -A;
        dxdt.I() = A - B;
        dxdt.R() = B;
    }
};

}  // namespace sir_intexp_reparam
}  // namespace country
}  // namespace epidemics
