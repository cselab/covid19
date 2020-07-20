#include "base.h"

namespace epidemics {
namespace cantons {
namespace seii_c {

/// State of the system.
template <typename T>
struct State : StateBase<T, 4> {
    using StateBase<T, 4>::StateBase;

    T &S (size_t i) { return this->v_[0 * this->numRegions_ + i]; }
    T &E (size_t i) { return this->v_[1 * this->numRegions_ + i]; }
    T &Ir(size_t i) { return this->v_[2 * this->numRegions_ + i]; }
    T &Iu(size_t i) { return this->v_[3 * this->numRegions_ + i]; }

    const T &S (size_t i) const { return this->v_[0 * this->numRegions_ + i]; }
    const T &E (size_t i) const { return this->v_[1 * this->numRegions_ + i]; }
    const T &Ir(size_t i) const { return this->v_[2 * this->numRegions_ + i]; }
    const T &Iu(size_t i) const { return this->v_[3 * this->numRegions_ + i]; }
};

template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 5;

    T beta;   /// Transmission rate.
    T nu;     /// Corrective multiplicative factor for Cij.
    T alpha;  /// Fraction of documented infections.
    T Z;      /// Average latency period.
    T D;      /// Average duration of infection.
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
        const double * __restrict__ invNi = dp_.invNi.data();
        for (size_t i = 0; i < dp_.numRegions; ++i) {
            T sumIC_N = ZERO + dp_.getExternalCommutersIu(day, i);
            for (size_t j = 0; j < dp_.numRegions; ++j)
                sumIC_N += x.Iu(j) * this->C_plus_Ct(i, j) * invNi[j];
            // printf("i=%zu invNi=%lg sumIC_N=%lg sum_SC_N=%lg\n", i, invNi[i], sumIC_N, sum_SC_N);
            T A = p.beta * x.S(i) * invNi[i] * (x.Iu(i) + p.nu * sumIC_N);
            T E_Z = x.E(i) / p.Z;

            T dS = -A;
            T dE = A - E_Z;
            T dIr = p.alpha * E_Z - x.Ir(i) / p.D;
            T dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;

            dxdt.S(i) = dS;
            dxdt.E(i) = dE;
            dxdt.Ir(i) = dIr;
            dxdt.Iu(i) = dIu;
        }
    }
};

}  // namesapce seii_c
}  // namespace cantons
}  // namespace epidemics
