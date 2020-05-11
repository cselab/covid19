#include "common.h"

namespace epidemics {
namespace cantons {
namespace seii_c {

/// Human-friendly wrapper around RawState.
struct State : StateBase<4> {
    using StateBase<4>::StateBase;

    double &S (size_t i) { return v_[0 * numRegions_ + i]; }
    double &E (size_t i) { return v_[1 * numRegions_ + i]; }
    double &Ir(size_t i) { return v_[2 * numRegions_ + i]; }
    double &Iu(size_t i) { return v_[3 * numRegions_ + i]; }

    double S (size_t i) const { return v_[0 * numRegions_ + i]; }
    double E (size_t i) const { return v_[1 * numRegions_ + i]; }
    double Ir(size_t i) const { return v_[2 * numRegions_ + i]; }
    double Iu(size_t i) const { return v_[3 * numRegions_ + i]; }
};

struct Parameters {
    double beta;   // Transmission rate.
    double nu;     // Corrective multiplicative factor for Cij.
    double alpha;  // Fraction of documented infections.
    double Z;      // Average latency period.
    double D;      // Average duration of infection.
};

struct Solver : SolverBase<Solver, State, Parameters> {
    using SolverBase<Solver, State, Parameters>::SolverBase;

    void rhs(int day, Parameters p, const State &x, State &dxdt) const;
};

}  // namesapce seii_c
}  // namespace cantons
}  // namespace epidemics
