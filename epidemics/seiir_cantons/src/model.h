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

struct MultiSEIIRState {
    // Per-region variables.
    std::vector<double> S;
    std::vector<double> E;
    std::vector<double> Idoc;
    std::vector<double> Iundoc;
    std::vector<double> N;

    void reset();
};

std::vector<MultiSEIIRState> solve_seiir(
        Parameters parameters, MultiSEIIRState initialState, int days);


class RegionalSEIIR {
public:
    RegionalSEIIR(std::vector<int> population, std::vector<int> commuteMatrix);

    double M(int from, int to) const {
        return M_[from * numRegions_ + to];
    }

    void deterministicRHS(const Parameters &p, const MultiSEIIRState &x, MultiSEIIRState &dxdt) const;

    std::vector<MultiSEIIRState> solve(Parameters parameters, int days) const;
private:

    size_t numRegions_;
    std::vector<int> N_;  // Population.
    std::vector<int> M_;  // Flattened matrix Mij.
    int R_;               // Number of regions.
};
