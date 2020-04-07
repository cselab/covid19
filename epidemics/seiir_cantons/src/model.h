#pragma once

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
 * Using custom types with boost::odeint is not that simple. Instead, we use a
 * single std::vector<double> to store the whole state of the simulation.
 */
using State = std::vector<double>;

/*
 * The view classes for reading and writing the state vector in a human-readable way.
 */
struct MultiRegionStateView {
    double &S (size_t i) const { return p_[0 * numRegions_ + i]; }
    double &E (size_t i) const { return p_[1 * numRegions_ + i]; }
    double &Ir(size_t i) const { return p_[2 * numRegions_ + i]; }
    double &Iu(size_t i) const { return p_[3 * numRegions_ + i]; }
    double &N (size_t i) const { return p_[4 * numRegions_ + i]; }

    MultiRegionStateView(size_t numRegions, double *p) :
        numRegions_{numRegions}, p_{p}
    { }
    MultiRegionStateView(State &state) :
        MultiRegionStateView{state.size() / 5, state.data()}
    { }

private:
    size_t numRegions_;
    double *p_;
};

struct MultiRegionStateConstView {
    double S (size_t i) const { return p_[0 * numRegions_ + i]; }
    double E (size_t i) const { return p_[1 * numRegions_ + i]; }
    double Ir(size_t i) const { return p_[2 * numRegions_ + i]; }
    double Iu(size_t i) const { return p_[3 * numRegions_ + i]; }
    double N (size_t i) const { return p_[4 * numRegions_ + i]; }

    MultiRegionStateConstView(size_t numRegions, const double *p) :
        numRegions_{numRegions}, p_{p}
    {}
    MultiRegionStateConstView(const State &state) :
        MultiRegionStateConstView{state.size() / 5, state.data()}
    { }

private:
    size_t numRegions_;
    const double *p_;
};

State createEmptyState(int numRegions);

class Solver {
public:
    Solver(std::vector<double> commuteMatrix);

    double M(int from, int to) const {
        return M_[from * numRegions_ + to];
    }
    double Mt(int to, int from) const {
        return Mt_[to * numRegions_ + from];
    }

    std::vector<State> solve(
            Parameters parameters,
            State initialState,
            int days) const;

private:
    void deterministicRHS(Parameters p, const State &x, State &dxdt) const;

    size_t numRegions_;
    std::vector<double> M_;   // Flattened matrix Mij.
    std::vector<double> Mt_;  // Flattened matrix Mji.
    int R_;                  // Number of regions.
};
