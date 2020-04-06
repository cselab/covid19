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
using MultiSEIIRState = std::vector<double>;

std::vector<MultiSEIIRState> solve_seiir(
        Parameters parameters, MultiSEIIRState initialState, int days);


class MultiSEIIR {
public:
    MultiSEIIR(std::vector<double> commuteMatrix);

    double M(int from, int to) const {
        return M_[from * numRegions_ + to];
    }
    double Mt(int to, int from) const {
        return Mt_[to * numRegions_ + from];
    }

    std::vector<MultiSEIIRState> solve(
            Parameters parameters,
            MultiSEIIRState initialState,
            int days) const;

private:
    void deterministicRHS(Parameters p, const MultiSEIIRState &x, MultiSEIIRState &dxdt) const;

    size_t numRegions_;
    std::vector<double> M_;   // Flattened matrix Mij.
    std::vector<double> Mt_;  // Flattened matrix Mji.
    int R_;                  // Number of regions.
};
