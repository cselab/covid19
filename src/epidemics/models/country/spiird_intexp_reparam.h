#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace spiird_intexp_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 7;

    T R0;     /// Reproduction number.
    T D;      /// Recovery time.
    T Y;      /// Preasymptomatic time.
    T alpha;  /// Reported rate
    T eps;    /// Mortality rate
    T tact;   /// Day of intervention.
    T k;      /// Exponential decay.
};

/// SEIIR state has 6 elements: S, P, Ir, Iu, R, D.
template <typename T>
struct State : StateBase<T, 6> {
    using StateBase<T, 6>::StateBase;

    T &S()  { return this->v_[0]; }
    T &P()  { return this->v_[1]; }
    T &Ir() { return this->v_[2]; }
    T &Iu() { return this->v_[3]; }
    T &R()  { return this->v_[4]; }
    T &D()  { return this->v_[5]; }

    const T &S()  const { return this->v_[0]; }
    const T &P()  const { return this->v_[1]; }
    const T &Ir() const { return this->v_[2]; }
    const T &Iu() const { return this->v_[3]; }
    const T &R()  const { return this->v_[4]; }
    const T &D()  const { return this->v_[5]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {


        auto invY   = 1 / p.Y;
        auto invD   = 1 / p.D;
        double invN = 1. / dp_.N;

        T r0 = intervention_exp(p.R0, t, p.k, p.tact);

        auto C1 = invN * r0 * invD * x.S() * (x.Iu() + x.P());
        auto C2 = invY * x.P();
        auto C3 = invD * x.Ir();
        auto C4 = invD * x.Iu();
            
        dxdt.S()  = -C1;
        dxdt.P()  = C1 - C2;
        dxdt.Ir() = p.alpha * C2 - C3;
        dxdt.Iu() = (1.0-p.alpha) * C2 - C4;
        dxdt.R()  = (1.0-p.eps) * C3 + C4;
        dxdt.D()  = p.eps * C3;
    }
};

}  // namespace spiird_intexp_reparam
}  // namespace country
}  // namespace epidemics
