#include "model.h"
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using Stepper = boost::numeric::odeint::runge_kutta_dopri5<State>;


static std::vector<double> transposeMatrix(const std::vector<double> &m, size_t N) {
    std::vector<double> out(N * N, 0.0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            out[i * N + j] = m[j * N + i];
    return out;
}

State createEmptyState(int numRegions) {
    return State(5 * numRegions, 0.0);
}

Solver::Solver(std::vector<double> commuteMatrix) :
    numRegions_{static_cast<size_t>(std::sqrt(commuteMatrix.size()) + 0.1)},
    M_{std::move(commuteMatrix)},
    Mt_{transposeMatrix(M_, numRegions_)}
{
    if (M_.size() != numRegions_ * numRegions_) {
        fprintf(stderr, "The number of elements of Mij is not a perfect square: %zu\n", M_.size());
        exit(1);
    }
}

std::vector<State> Solver::solve(Parameters parameters, State initialState, int days) const {
    const size_t NUM_VARS = 5;
    if (initialState.size() != numRegions_ * NUM_VARS) {
        fprintf(stderr, "Expected %zu elements in initialState, got %zu\n",
            numRegions_ * NUM_VARS, initialState.size());
        exit(1);
    }

    const int STEPS_PER_DAY = 10;
    const double dt = 1.0 / STEPS_PER_DAY;

    auto rhs = [this, parameters](const State &x, State &dxdt, double /*t*/) {
        deterministicRHS(parameters, x, dxdt);
    };

    std::vector<State> result;

    // Observer gets called for each time step evaluated by the integrator.
    // We consider only every `STEPS_PER_DAY` steps (skipping the first one as well).
    auto observer = [&result, cnt = 0](const State &y, double /*t*/) mutable {
        if (cnt % STEPS_PER_DAY == 0 && cnt > 0)
            result.push_back(y);
        ++cnt;
    };

    boost::numeric::odeint::integrate_n_steps(
            Stepper{}, rhs, initialState, 0.0, dt, days * STEPS_PER_DAY, observer);
    return result;
}

void Solver::deterministicRHS(
        Parameters p,
        const State &x_,
        State &dxdt_) const {
    if (dxdt_.size() != x_.size())
        throw std::runtime_error("dxdt does not have the expected size.");
    if (x_.size() != numRegions_ * 5)
        throw std::runtime_error("x does not have the expected size.");

    MultiRegionStateConstView x{numRegions_, x_.data()};
    MultiRegionStateView dxdt{numRegions_, dxdt_.data()};

    for (size_t i = 0; i < numRegions_; ++i) {
        double A = p.beta * x.S(i) / x.N(i) * x.Ir(i);
        double B = p.beta * x.S(i) / x.N(i) * p.mu * x.Iu(i);
        double E_Z = x.E(i) / p.Z;

        double dS = -(A + B);
        double dE = A + B - E_Z;
        double dIr = p.alpha * E_Z - x.Ir(i) / p.D;
        double dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;
        double dN = 0;

        double inv = 1 / (x.N(i) - x.Ir(i));
        for (size_t j = 0; j < numRegions_; ++j) {
            double Tij = this->M (i, j) / (x.N(j) - x.Ir(j));
            double Tji = this->Mt(i, j) * inv;
            dS += p.theta * (Tij * x.S(j) - Tji * x.S(i));
            dE += p.theta * (Tij * x.E(j) - Tji * x.E(i));
            // Documented infected people are in quarantine, they do not move around.
            // dIr += p.theta * (Tij * x.Ir(j) - Tji * x.Ir(i));
            dIu += p.theta * (Tij * x.Iu(j) - Tji * x.Iu(i));
            dN += p.theta * (this->M(i, j) - this->Mt(i, j));
        }

        dxdt.S(i) = dS;
        dxdt.E(i) = dE;
        dxdt.Ir(i) = dIr;
        dxdt.Iu(i) = dIu;
        dxdt.N(i) = dN;
    }
}
