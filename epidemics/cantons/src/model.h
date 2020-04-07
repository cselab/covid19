#pragma once

#include <cassert>
#include <vector>

/// Parameters to be optimized for.
struct Parameters {
    double beta;   // Transmission rate.
    double mu;     // Reduction factor for transmission rate of undocumented individuals.
    double alpha;  // Fraction of documented infections.
    double Z;      // Average latency period.
    double D;      // Average duration of infection.
    double theta;  // Corrective multiplicative factor for Mij.

    double computeRe() const {
        return alpha * beta * D + (1 - alpha) * mu * beta * D;
    }
};

/*
 * Using custom types with boost::odeint is not that simple.
 * Instead, we use a single std::vector<double> to store the
 * whole state of the simulation, and the wrapper below to
 * access the data in a human-friendly way.
 */
using RawState = std::vector<double>;

/// Human-friendly wrapper around RawState.
struct State {
    static constexpr int kVarsPerRegion = 5;

    /// Create a state with all values set to 0.
    explicit State(size_t numRegions) :
        numRegions_{numRegions},
        v_(kVarsPerRegion * numRegions, 0.0)
    { }
    explicit State(RawState state) :
        numRegions_{state.size() / kVarsPerRegion},
        v_{std::move(state)}
    {
        assert(v_.size() % kVarsPerRegion == 0);
    }
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

    size_t numRegions() const noexcept { return numRegions_; }
    const RawState &raw() const & { return v_; }
    RawState raw() && { return std::move(v_); }

private:
    size_t numRegions_;
    RawState v_;
};

class Solver {
public:
    Solver(std::vector<double> commuteMatrix);

    double M(int from, int to) const {
        return M_[from * numRegions_ + to];
    }
    double Mt(int to, int from) const {
        return Mt_[to * numRegions_ + from];
    }

    std::vector<State> solve(const Parameters &parameters, State initialState, int days) const;
    std::vector<State> solve(const Parameters &parameters, RawState initialState, int days) const;

private:
    void deterministicRHS(Parameters p, const RawState &x, RawState &dxdt) const;

    size_t numRegions_;
    std::vector<double> M_;   // Flattened matrix Mij.
    std::vector<double> Mt_;  // Flattened matrix Mji.
    int R_;                  // Number of regions.
};
