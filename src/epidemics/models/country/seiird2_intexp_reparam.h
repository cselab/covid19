#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace seiird2_intexp_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 8;

    T R0;     /// Reproduction number.
    T D;      /// Recovery time.
    T Z;      /// Incubation time.
    T mu;     /// Transmissivity modulation factor
    T alpha;  /// Reported rate
    T eps;    /// Mortality rate
    T tact;   /// Day of intervention.
    T k;      /// Expoenntial decay.
};

/// SEIIR state has 6 elements: S, E, Ir, Iu, R, D.
template <typename T>
struct State : StateBase<T, 6> {
    using StateBase<T, 6>::StateBase;

    T &S()  { return this->v_[0]; }
    T &E()  { return this->v_[1]; }
    T &Ir() { return this->v_[2]; }
    T &Iu() { return this->v_[3]; }
    T &R()  { return this->v_[4]; }
    T &D()  { return this->v_[5]; }

    const T &S()  const { return this->v_[0]; }
    const T &E()  const { return this->v_[1]; }
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


        auto invZ   = 1 / p.Z;
        auto invD   = 1 / p.D;
        double invN = 1. / dp_.N;

        T r0 = intervention_exp(p.R0, t, p.k, p.tact);

        auto C1 = invN * r0 * invD * x.S() * x.Ir();
        auto C2 = invN * r0 * invD * x.S() * (p.mu * x.Iu());
        auto C3 = p.alpha * (invZ * x.E());
        auto C4 = (1 - p.alpha) * (invZ * x.E());
        auto C5 = invD * x.Ir();
        auto C6 = invD * x.Iu();
            
        dxdt.S()  = -(C1 + C2);
        dxdt.E()  = C1 + C2 - C3 - C4;
        dxdt.Ir() = C3 - C5;
        dxdt.Iu() = C4 - C6;
        dxdt.R()  = (1.0-p.eps)*C5 + C6;
        dxdt.D()  = p.eps*C5;
    }
};

}  // namespace seiird2_intexp_reparam
}  // namespace country
}
