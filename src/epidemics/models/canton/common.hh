#pragma once

#include "common.h"
#include <epidemics/utils/signal.h>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cstdio>

using Stepper = boost::numeric::odeint::runge_kutta_dopri5<RawState>;

template <typename Derived, typename State, typename Parameters>
std::vector<State>
SolverBase<Derived, State, Parameters>::solve(
        const Parameters &parameters,
        RawState initialState,
        int days) const
{
    if (initialState.size() != modelData_.numRegions * State::kVarsPerRegion) {
        fprintf(stderr,
                "Expected %zu elements in initialState, got %zu\n",
                modelData_.numRegions * State::kVarsPerRegion,
                initialState.size());
        exit(1);
    }

    const int STEPS_PER_DAY = 10;
    const double dt = 1.0 / STEPS_PER_DAY;
    std::vector<State> result;

    // Observer gets called for each time step evaluated by the integrator.
    // We consider only every `STEPS_PER_DAY` steps (skipping the first one as well).
    auto observer = [&result, cnt = 0, verbose = this->verbose_](
            const RawState &y,
            double /*t*/) mutable {
        if (check_signals_func)
            check_signals_func();
        if (cnt % STEPS_PER_DAY == 0 && cnt > 0) {
            result.push_back(State{y});
            if (verbose) {
                double total_Ir = 0;
                for (size_t i = 0; i < result.back().numRegions(); ++i)
                    total_Ir += result.back().Ir(i);
                printf("Total Ir on day %zu: %lg\n", result.size(), total_Ir);
            }
        }
        ++cnt;
    };

    auto rhs = [this, parameters](const RawState &x_, RawState &dxdt_, double t) {
        int day = static_cast<int>(t);
        if (dxdt_.size() != x_.size())
            throw std::runtime_error("dxdt does not have the expected size.");
        if (x_.size() != modelData_.numRegions * State::kVarsPerRegion)
            throw std::runtime_error("x does not have the expected size.");

        // This is a tricky part, we transform RawState to State during
        // computation and then at the end transform it back.
        State x{std::move(const_cast<RawState&>(x_))};
        State dxdt{std::move(dxdt_)};

        derived()->rhs(day, parameters, x, dxdt);

        /// Transform State back to RawState.
        const_cast<RawState &>(x_) = std::move(x).raw();
        dxdt_ = std::move(dxdt).raw();
    };

    boost::numeric::odeint::integrate_n_steps(
            Stepper{}, rhs, initialState, 0.0, dt, days * STEPS_PER_DAY, observer);
    return result;
}
