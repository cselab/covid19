#pragma once

#include "base.h"

namespace epidemics {
namespace country {
namespace seiir_int {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 8;

    T beta;   /// Transmission rate.
    T mu;     /// Proportionality parameter.
    T alpha;  /// Fraction of reported incidents.
    T Z;      /// Average latency period.
    T D;      /// Average duration of infection.
    T tact;   /// Day of intervention.
    T dtact;  /// Duration of intervention.
    T kbeta;  /// Multiplicator beta after intervention.
};

/// SEIIR state has 5 elements: S, E, Ir, Iu, R.
template <typename T>
struct State : StateBase<T, 5> {
    using StateBase<T, 5>::StateBase;

    T &S()  { return this->v_[0]; }
    T &E()  { return this->v_[1]; }
    T &Ir() { return this->v_[2]; }
    T &Iu() { return this->v_[3]; }
    T &R()  { return this->v_[4]; }

    const T &S()  const { return this->v_[0]; }
    const T &E()  const { return this->v_[1]; }
    const T &Ir() const { return this->v_[2]; }
    const T &Iu() const { return this->v_[3]; }
    const T &R()  const { return this->v_[4]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {

        T beta;
        if (t < p.tact - 0.5*p.dtact) {
           beta = p.beta;
        } else if (t < p.tact + 0.5*p.dtact) {
           beta = (1. - (t - 0.5*p.dtact - p.tact) / p.dtact * (1. - p.kbeta)) * p.beta;
        } else {
           beta = p.kbeta * p.beta;
        }

        auto invZ   = 1 / p.Z;
        auto invD   = 1 / p.D;
        double invN = 1. / data_.N;

        auto C1 = invN * beta * x.S() * x.Ir();
        auto C2 = invN * beta * x.S() * (p.mu * x.Iu());
        auto C3 = p.alpha * (invZ * x.E());
        auto C4 = (1 - p.alpha) * (invZ * x.E());

        dxdt.S()  = -(C1 + C2);
        dxdt.E()  = C1 + C2 - invZ * x.E();
        dxdt.Ir() = C3 - invD * x.Ir();
        dxdt.Iu() = C4 - invD * x.Iu();
        dxdt.R()  = -dxdt.S() - dxdt.E() - dxdt.Ir() - dxdt.Iu();
    }
};

}  // namespace seiir_int
}  // namespace country
}  // namespace epidemics
