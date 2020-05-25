#pragma once

#include "base.h"

namespace epidemics {
namespace cantons {
namespace seiin {

/// State of the system.
template <typename T>
struct State : StateBase<T, 5> {
    using StateBase<T, 5>::StateBase;
    T &S (size_t i) { return this->v_[0 * this->numRegions_ + i]; }
    T &E (size_t i) { return this->v_[1 * this->numRegions_ + i]; }
    T &Ir(size_t i) { return this->v_[2 * this->numRegions_ + i]; }
    T &Iu(size_t i) { return this->v_[3 * this->numRegions_ + i]; }
    T &N (size_t i) { return this->v_[4 * this->numRegions_ + i]; }

    const T &S (size_t i) const { return this->v_[0 * this->numRegions_ + i]; }
    const T &E (size_t i) const { return this->v_[1 * this->numRegions_ + i]; }
    const T &Ir(size_t i) const { return this->v_[2 * this->numRegions_ + i]; }
    const T &Iu(size_t i) const { return this->v_[3 * this->numRegions_ + i]; }
    const T &N (size_t i) const { return this->v_[4 * this->numRegions_ + i]; }
};


/// Lightweight parameters (optimized for).
template <typename T>
struct Parameters {
    static constexpr size_t numParameters = 6;

    T beta;   /// Transmission rate.
    T mu;     /// Reduction factor for transmission rate of undocumented individuals.
    T alpha;  /// Fraction of documented infections.
    T Z;      /// Average latency period.
    T D;      /// Average duration of infection.
    T theta;  /// Corrective multiplicative factor for Mij.
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
        for (size_t i = 0; i < modelData_.numRegions; ++i) {
            double extComIu = modelData_.getExternalCommutersIu(day, i);
            T A = p.beta * x.S(i) / x.N(i) * (x.Ir(i) + extComIu);
            T B = p.beta * x.S(i) / x.N(i) * p.mu * x.Iu(i);
            T E_Z = x.E(i) / p.Z;

            T dS = -(A + B);
            T dE = A + B - E_Z;
            T dIr = p.alpha * E_Z - x.Ir(i) / p.D;
            T dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;
            T dN = ZERO;

            T inv = 1 / (x.N(i) - x.Ir(i));
            // NOTE: nonzero_Mij will not work if non-zero elements do not on symmetric places.

            //for (size_t j = 0; j < modelData_.numRegions; ++j) { // XXX
            for (size_t j : this->nonzero_Mij(i)) {
                auto Tij = this->M(i, j) / (x.N(j) - x.Ir(j));
                auto Tji = this->M(j, i) * inv;
                dS += p.theta * (Tij * x.S(j) - Tji * x.S(i));
                dE += p.theta * (Tij * x.E(j) - Tji * x.E(i));
                // Documented infected people are in quarantine, they do not move around.
                // dIr += p.theta * (Tij * x.Ir(j) - Tji * x.Ir(i));
                dIu += p.theta * (Tij * x.Iu(j) - Tji * x.Iu(i));
                dN += p.theta * (this->M(i, j) - this->M(j, i));
            }

            dxdt.S(i) = dS;
            dxdt.E(i) = dE;
            dxdt.Ir(i) = dIr;
            dxdt.Iu(i) = dIu;
            dxdt.N(i) = dN;
        }
    }
};

}  // namespace seiin
}  // namespace cantons
}  // namespace epidemics
