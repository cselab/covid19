#include "common.h"

namespace sei_c {

/// Human-friendly wrapper around RawState.
struct State : StateBase<3> {
    using StateBase<3>::StateBase;

    double &S(size_t i) { return v_[0 * numRegions_ + i]; }
    double &E(size_t i) { return v_[1 * numRegions_ + i]; }
    double &I(size_t i) { return v_[2 * numRegions_ + i]; }

    double S(size_t i) const { return v_[0 * numRegions_ + i]; }
    double E(size_t i) const { return v_[1 * numRegions_ + i]; }
    double I(size_t i) const { return v_[2 * numRegions_ + i]; }
    // XXX needed by `observer` in common.hh
    double Ir(size_t i) const { return I(i); }
};

struct Parameters {
    double beta;   // Transmission rate.
    double nu;     // Corrective multiplicative factor for Cij.
    double Z;      // Average latency period.
    double D;      // Average duration of infection.
    double tact;   // Time of intervention.
    double kbeta;   // Factor for `beta` after intervention.
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    void rhs(int day, Parameters p, const State &x, State &dxdt) const;
};

}  // namesapce sei_c
