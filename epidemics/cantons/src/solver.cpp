#include "solver.h"
#include "utils.h"
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using Stepper = boost::numeric::odeint::runge_kutta_dopri5<RawState>;


static std::vector<double> transposeMatrix(const std::vector<double> &m, size_t N) {
    std::vector<double> out(N * N, 0.0);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            out[i * N + j] = m[j * N + i];
    return out;
}

Solver::Solver(ModelData modelData) :
    modelData_{std::move(modelData)}
{ }

std::vector<State> Solver::solve(const Parameters &parameters, State initialState, int days) const {
    return solve(parameters, std::move(initialState).raw(), days);
}

std::vector<State> Solver::solve(const Parameters &parameters, RawState initialState, int days) const {
    if (initialState.size() != modelData_.numRegions * State::kVarsPerRegion) {
        fprintf(stderr, "Expected %zu elements in initialState, got %zu\n",
            modelData_.numRegions * State::kVarsPerRegion, initialState.size());
        exit(1);
    }

    const int STEPS_PER_DAY = 10;
    const double dt = 1.0 / STEPS_PER_DAY;
    std::vector<State> result;

    // Observer gets called for each time step evaluated by the integrator.
    // We consider only every `STEPS_PER_DAY` steps (skipping the first one as well).
    auto observer = [&result, cnt = 0](const RawState &y, double /*t*/) mutable {
        if (cnt % STEPS_PER_DAY == 0 && cnt > 0)
            result.push_back(State{y});
        ++cnt;
    };

    auto rhs = [this, parameters](const RawState &x, RawState &dxdt, double t) {
        int day = static_cast<int>(t);
        // const_cast is (kind of) valid, because we move from and move back to
        // x, so no change is made in the end.
        deterministicRHS(day, parameters, const_cast<RawState &>(x), dxdt);
    };

    boost::numeric::odeint::integrate_n_steps(
            Stepper{}, rhs, initialState, 0.0, dt, days * STEPS_PER_DAY, observer);
    return result;
}

void Solver::deterministicRHS(
        int day,
        Parameters p,
        const RawState &x_,
        RawState &dxdt_) const {
    if (dxdt_.size() != x_.size())
        throw std::runtime_error("dxdt does not have the expected size.");
    if (x_.size() != modelData_.numRegions * 5)
        throw std::runtime_error("x does not have the expected size.");

    // This is a tricky part, we transform RawState to State during
    // computation and then at the end transform it back.
    State x{std::move(const_cast<RawState&>(x_))};
    State dxdt{std::move(dxdt_)};

    for (size_t i = 0; i < modelData_.numRegions; ++i) {
        double external_cases = modelData_.getExternalCasesAt(day, i);
        double A = p.beta * x.S(i) / x.N(i) * (x.Ir(i) + external_cases);
        double B = p.beta * x.S(i) / x.N(i) * p.mu * x.Iu(i);
        double E_Z = x.E(i) / p.Z;

        double dS = -(A + B);
        double dE = A + B - E_Z;
        double dIr = p.alpha * E_Z - x.Ir(i) / p.D;
        double dIu = E_Z - p.alpha * E_Z - x.Iu(i) / p.D;
        double dN = 0;

        double inv = 1 / (x.N(i) - x.Ir(i));
        for (size_t j = 0; j < modelData_.numRegions; ++j) {
            double Tij = this->M(i, j) / (x.N(j) - x.Ir(j));
            double Tji = this->M(j, i) * inv;
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

    /// Transform State back to RawState.
    const_cast<RawState &>(x_) = std::move(x).raw();
    dxdt_ = std::move(dxdt).raw();
}
