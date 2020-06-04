#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace seir {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 3;

    T beta;   /// Transmission rate.
    T gamma;  /// Recovery rate.
    T a;      /// Inverse of average incubation period.
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
    void rhs(double /* t */,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {
        double invN = 1. / data_.N;
        auto A = invN * p.beta * x.I() * x.S();
        auto B = p.a*x.E();
        auto C = p.gamma * x.I();

        dxdt.S() = -A;
        dxdt.E() = +A-B;
        dxdt.I() = +B-C;
        dxdt.R() = +C;
    }
};

}  // namespace seir
}  // namespace country
}  // namespace epidemics
