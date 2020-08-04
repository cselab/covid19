#pragma once

#include "base.h"
#include "intervention.h"

namespace epidemics {
namespace country {
namespace cz_intstep {

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 14;

    T R0;     /// Reproduction number.
    T gamma;  /// Recovery rate.
    T sigma;  /// Latency rate.
    T eps1;   ///
    T eps2;   ///
    T eps3;   ///
    T eps4;   ///
    T omega1; ///
    T omega2; ///
    T omega3; ///
    T omega4; ///
    T omega5; ///
    T tact;   /// Day of intervention.
    T kbeta;  /// Multiplicator.
};

/// CZ has 10 elements
template <typename T>
struct State : StateBase<T, 10> {
    using StateBase<T, 10>::StateBase;

    T &S()  { return this->v_[0]; }
    T &E()  { return this->v_[1]; }
    T &I()  { return this->v_[2]; }
    T &P()  { return this->v_[3]; }
    T &H1() { return this->v_[4]; }
    T &H2() { return this->v_[5]; }
    T &U()  { return this->v_[6]; }
    T &R()  { return this->v_[7]; }
    T &D()  { return this->v_[8]; }
    T &C()  { return this->v_[9]; }

    const T &S()  const { return this->v_[0]; }
    const T &E()  const { return this->v_[1]; }
    const T &I()  const { return this->v_[2]; }
    const T &P()  const { return this->v_[3]; }
    const T &H1() const { return this->v_[4]; }
    const T &H2() const { return this->v_[5]; }
    const T &U()  const { return this->v_[6]; }
    const T &R()  const { return this->v_[7]; }
    const T &D()  const { return this->v_[8]; }
    const T &C()  const { return this->v_[9]; }
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> &x,
             State<T> & __restrict__ dxdt) const
    {

        T r0 = intervention_step(p.R0, t, p.kbeta, p.tact);

        double beta = r0 / dp_.N * p.gamma;

        auto A = beta * x.I() * x.S();
        auto B = p.sigma * x.E();
        auto C = p.gamma * x.I();
        auto D = p.omega1 * x.P();
        auto E = p.omega2 * x.H1();
        auto F = p.omega3 * x.H2();
        auto G = p.omega4 * x.U();
        auto H = p.omega5 * x.U();

        dxdt.S()  = -A;
        dxdt.E()  = +A-B;
        dxdt.I()  = +B-C;
        dxdt.P()  = +p.eps1*C-D;
        dxdt.H1() = +D-E;
        dxdt.H2() = +(1.0-p.eps2)*E-F;
        dxdt.U()  = +p.eps2*E-(1.0-p.eps4)*G-p.eps4*H;
        dxdt.R()  = +(1.0-p.eps1)*C+(1.0-p.eps3)*F+(1.0-p.eps4)*G;
        dxdt.D()  = +p.eps2*F+p.eps4*H;
        dxdt.C()  = +C;
    }
};

}  // namespace cz_intstep
}  // namespace country
}  // namespace epidemics
