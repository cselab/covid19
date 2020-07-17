#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace seir_intexp_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 5;

    T R0;     /// Reproduction number.
    T D;      /// Recovery time.
    T Z;      /// Incubation period.
    T tact;   /// Day of intervention.
    T k;      /// Exponential decay.
};

/// SIR state has 3 elements: S, E, I, R.
template <typename T>
struct State : StateBase<T, 4> {
    using StateBase<T, 4>::StateBase;

    T &S() { return this->v_[0]; }
    T &E() { return this->v_[1]; }
    T &I() { return this->v_[2]; }
    T &R() { return this->v_[3]; }

    const T &S() const { return this->v_[0]; }
    const T &E() const { return this->v_[1]; }
    const T &I() const { return this->v_[2]; }
    const T &R() const { return this->v_[3]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {

        T r0;
        if (t <= p.tact) {
           r0 = p.R0;
        } else {
           using std::exp;
           double dt = p.tact-t; // careful: diff will not work
           r0 = p.R0*exp(p.k*dt);
        }

        double invN = 1. / dp_.N;
        auto invD = 1. / p.D;
        auto invZ = 1. / p.Z;

        auto A = invN * r0 * invD * x.I() * x.S();
        auto B = invZ * x.E();
        auto C = invD * x.I();

        dxdt.S() = -A;
        dxdt.E() = +A-B;
        dxdt.I() = +B-C;
        dxdt.R() = +C;
    }
};

}  // namespace seir_int
}  // namespace country
}  // namespace epidemics
