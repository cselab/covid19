#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace sir_intrem_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 4;

    T R0;     /// Transmission rate.
    T D;      /// Recovery time.
    T tact;   /// Day of intervention.
    T l;      /// Removal fraction.
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
        auto invD = 1. / p.D;

        T removal;

        if(t > p.tact) {
            removal = p.l*x.S();
        }
        else
        {
            removal = 0.0*x.S();
        }

        auto A = invN * p.R0 * invD * x.I() * x.S();
        auto B = invD * x.I();

        dxdt.S() = -A - removal;
        dxdt.I() = A - B;
        dxdt.R() = B + removal;
    }
};

}  // namespace sir_intexp_reparam
}  // namespace country
}  // namespace epidemics
