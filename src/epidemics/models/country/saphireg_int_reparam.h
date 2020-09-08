#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace saphireg_int_reparam {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 11;

    T R0;     /// Reproduction number.
    T D;      /// Recovery time.
    T F;      /// Removal time.
    T Y;      /// Preasymptomatic time.
    T Z;      /// Latency period.
    T mu;     /// Transmissivity modulation factor
    T alpha;  /// Reported rate
    T eps;    /// Mortality rate
    T tact;   /// Day of intervention.
    T dtact;  /// Duration of intervention.
    T kbeta;  /// Multiplicator beta after intervention.
};

/// SAPHIRE state has 7 elements: S, E P, Ir, Iu, R, D.
template <typename T>
struct State : StateBase<T, 7> {
    using StateBase<T, 7>::StateBase;

    T &S()  { return this->v_[0]; }
    T &E()  { return this->v_[1]; }
    T &P()  { return this->v_[2]; }
    T &Ir() { return this->v_[3]; }
    T &Iu() { return this->v_[4]; }
    T &R()  { return this->v_[5]; }
    T &D()  { return this->v_[6]; }

    const T &S()  const { return this->v_[0]; }
    const T &E()  const { return this->v_[1]; }
    const T &P()  const { return this->v_[2]; }
    const T &Ir() const { return this->v_[3]; }
    const T &Iu() const { return this->v_[4]; }
    const T &R()  const { return this->v_[5]; }
    const T &D()  const { return this->v_[6]; }
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
        auto invY   = 1 / p.Y;
        auto invD   = 1 / p.D;
        auto invF   = 1 / p.F;
        double invN = 1. / dp_.N;

        T r0 = intervention(p.R0, t, p.kbeta, p.tact, p.dtact);

        auto C1 = invN * r0 * invD * x.S() * (p.mu * (x.Iu() + x.P()) + x.Ir());
        auto C2 = invZ * x.E();
        auto C3 = invY * x.P();
        auto C4 = (1.0 - p.eps) * invD * x.Ir();
        auto C5 = invD * x.Iu();
        auto C6 = p.eps * invF * x.Ir();
            
        dxdt.S()  = -C1;
        dxdt.E()  = C1 - C2;
        dxdt.P()  = C2 - C3;
        dxdt.Ir() = p.alpha * C3 - C4 - C6;
        dxdt.Iu() = (1.0-p.alpha) * C3 - C5;
        dxdt.R()  = C4 + C5;
        dxdt.D()  = C6;
    }
};

}  // namespace saphireg_int_reparam
}  // namespace country
}  // namespace epidemics
