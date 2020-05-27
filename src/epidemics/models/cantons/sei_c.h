#include "base.h"

namespace epidemics {
namespace cantons {
namespace sei_c {

/// State of the system.
template <typename T>
struct State : StateBase<T, 3> {
    using StateBase<T, 3>::StateBase;

    T &S(size_t i) { return this->v_[0 * this->numRegions_ + i]; }
    T &E(size_t i) { return this->v_[1 * this->numRegions_ + i]; }
    T &I(size_t i) { return this->v_[2 * this->numRegions_ + i]; }

    const T &S(size_t i) const { return this->v_[0 * this->numRegions_ + i]; }
    const T &E(size_t i) const { return this->v_[1 * this->numRegions_ + i]; }
    const T &I(size_t i) const { return this->v_[2 * this->numRegions_ + i]; }
    // XXX needed by `observer` in base.hh
    const T &Ir(size_t i) const { return I(i); }
};

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 6;

    T beta;   /// Transmission rate.
    T nu;     /// Corrective multiplicative factor for Cij.
    T Z;      /// Average latency period.
    T D;      /// Average duration of infection.
    T tact;   /// Time of intervention.
    T kbeta;  /// Factor for `beta` after intervention.
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    template <typename T>
    void rhs(double t,
             Parameters<T> p,
             const State<T> & __restrict__ x,
             State<T> & __restrict__ dxdt) const
    {
        const T ZERO = 0 * p.beta;
        int day = static_cast<int>(t);
        T beta = p.beta;
        // FIXME: rhs should take `double day` (same in lambda `rhs`)
        //        then replace `day + 1` by `day`
        if (day + 1 > p.tact) {
          beta *= p.kbeta;
        }
        const double * __restrict__ invNi = modelData_.invNi.data();
        for (size_t i = 0; i < modelData_.numRegions; ++i) {
            T sumIC_N = ZERO;
            for (size_t j = 0; j < modelData_.numRegions; ++j) {
                sumIC_N += x.I(j) * this->C_plus_Ct(i, j) * invNi[j];
            }
            const double ext = modelData_.getExternalCommutersIu(day, i);
            const T beta_i = beta * (1 + modelData_.Ui[i]);
            const T A = beta_i * x.S(i) * invNi[i] * (
                    x.I(i) + p.nu * sumIC_N + ext);
            const T E_Z = x.E(i) / p.Z;

            const T dS = -A;
            const T dE = A - E_Z;
            const T dI = E_Z - x.I(i) / p.D;

            dxdt.S(i) = dS;
            dxdt.E(i) = dE;
            dxdt.I(i) = dI;
        }
    }
};

}  // namesapce sei_c
}  // namespace cantons
}  // namespace epidemics
