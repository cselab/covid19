#pragma once

#include "common.h"

namespace seiin {

/// Human-friendly wrapper around RawState.
struct State : StateBase<5> {
    using StateBase<5>::StateBase;
    double &S (size_t i) { return v_[0 * numRegions_ + i]; }
    double &E (size_t i) { return v_[1 * numRegions_ + i]; }
    double &Ir(size_t i) { return v_[2 * numRegions_ + i]; }
    double &Iu(size_t i) { return v_[3 * numRegions_ + i]; }
    double &N (size_t i) { return v_[4 * numRegions_ + i]; }

    double S (size_t i) const { return v_[0 * numRegions_ + i]; }
    double E (size_t i) const { return v_[1 * numRegions_ + i]; }
    double Ir(size_t i) const { return v_[2 * numRegions_ + i]; }
    double Iu(size_t i) const { return v_[3 * numRegions_ + i]; }
    double N (size_t i) const { return v_[4 * numRegions_ + i]; }
};


/// Lightweight parameters (optimized for).
struct Parameters {
    double beta;   // Transmission rate.
    double mu;     // Reduction factor for transmission rate of undocumented individuals.
    double alpha;  // Fraction of documented infections.
    double Z;      // Average latency period.
    double D;      // Average duration of infection.
    double theta;  // Corrective multiplicative factor for Mij.
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    void rhs(int day, Parameters p, const State &x, State &dxdt) const;
};

}  // namespace seiin
