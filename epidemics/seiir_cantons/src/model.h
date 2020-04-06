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
 *
 * The view classes below are used to read and write the state vector in a human-readable way.
 */
using MultiSEIIRState = std::vector<double>;

struct MultiSEIIRStateView {
    double &S (size_t i) const { return v_[0 * numRegions_ + i]; }
    double &E (size_t i) const { return v_[1 * numRegions_ + i]; }
    double &Ir(size_t i) const { return v_[2 * numRegions_ + i]; }
    double &Iu(size_t i) const { return v_[3 * numRegions_ + i]; }
    double &N (size_t i) const { return v_[4 * numRegions_ + i]; }

    MultiSEIIRStateView(std::vector<double> &v) : numRegions_{v.size() / 5}, v_{v} { }

private:
    size_t numRegions_;
    std::vector<double> &v_;
};

struct MultiSEIIRStateConstView {
    double S (size_t i) const { return v_[0 * numRegions_ + i]; }
    double E (size_t i) const { return v_[1 * numRegions_ + i]; }
    double Ir(size_t i) const { return v_[2 * numRegions_ + i]; }
    double Iu(size_t i) const { return v_[3 * numRegions_ + i]; }
    double N (size_t i) const { return v_[4 * numRegions_ + i]; }

    MultiSEIIRStateConstView(const std::vector<double> &v) : numRegions_{v.size() / 5}, v_{v} { }

private:
    size_t numRegions_;
    const std::vector<double> &v_;
};

std::vector<MultiSEIIRState> solve_seiir(
        Parameters parameters, MultiSEIIRState initialState, int days);


class MultiSEIIR {
public:
    MultiSEIIR(std::vector<int> population, std::vector<int> commuteMatrix);

    double M(int from, int to) const {
        return M_[from * numRegions_ + to];
    }

    std::vector<MultiSEIIRState> solve(
            Parameters parameters,
            MultiSEIIRState initialState,
            int days) const;

private:
    void deterministicRHS(const Parameters &p, const MultiSEIIRState &x, MultiSEIIRState &dxdt) const;

    size_t numRegions_;
    std::vector<int> N_;  // Population.
    std::vector<int> M_;  // Flattened matrix Mij.
    int R_;               // Number of regions.
};
