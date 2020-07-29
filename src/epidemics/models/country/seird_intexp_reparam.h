#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace seird_intexp_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 6;

    T R0;     /// Reproduction number.
    T D;      /// Recovery time.
    T Z;      /// Incubation period.
    T eps;    /// Death rate, from infected
    T tact;   /// Day of intervention.
    T k;      /// Exponential decay.
};

/// SEIRD state has 5 elements: S, E, I, R, D
template <typename T>
struct State : StateBase<T, 5> {
    using StateBase<T, 5>::StateBase;

    T &S() { return this->v_[0]; }
    T &E() { return this->v_[1]; }
    T &I() { return this->v_[2]; }
    T &R() { return this->v_[3]; }
    T &D() { return this->v_[4]; }

    const T &S() const { return this->v_[0]; }
    const T &E() const { return this->v_[1]; }
    const T &I() const { return this->v_[2]; }
    const T &R() const { return this->v_[3]; }
    const T &D() const { return this->v_[4]; }

};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {

        T r0 = intervention_exp(p.R0, t, p.k, p.tact);

        double invN = 1. / dp_.N;
        auto invD = 1. / p.D;
        auto invZ = 1. / p.Z;

        auto A = invN * r0 * invD * x.I() * x.S();
        auto B = invZ * x.E();
        auto C = invD * x.I();

        dxdt.S() = -A;
        dxdt.E() = +A-B;
        dxdt.I() = +B-C;
        dxdt.R() = (1.0-p.eps)*C;
        dxdt.D() = p.eps*C;

    }
};

}  // namespace seird_intexp
}  // namespace country
}  // namespace epidemics
