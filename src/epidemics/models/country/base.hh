#pragma once

#include "base.h"
#include <epidemics/utils/signal.h>

#include <boost/numeric/odeint.hpp>
#include <cstdio>

namespace epidemics {
namespace country {

template <typename Derived,
          template <typename> typename State,
          template <typename> typename Parameters>
template <typename T>
std::vector<State<T>>
SolverBase<Derived, State, Parameters>::solve(
        const Parameters<T> &parameters,
        typename State<T>::RawState initialState,
        const std::vector<double> &tEval) const
{
    using RawState = typename State<T>::RawState;
    using Stepper = boost::numeric::odeint::runge_kutta_dopri5<RawState>;

    const double dt = 1.0 / 10;
    std::vector<State<T>> result;
    result.reserve(tEval.size());

    auto observer = [&result](const RawState &y, double /*t*/) mutable {
        if (check_signals_func)
            check_signals_func();
        result.push_back(State<T>{y});
    };

    auto rhs = [this, parameters](const RawState &x_, RawState &dxdt_, double t) {
        // This is a tricky part, we transform RawState to State during
        // computation and then at the end transform it back.
        State<T> x{std::move(const_cast<RawState&>(x_))};
        State<T> dxdt{std::move(dxdt_)};

        derived()->rhs(t, parameters, x, dxdt);

        /// Transform State back to RawState.
        const_cast<RawState &>(x_) = std::move(x).raw();
        dxdt_ = std::move(dxdt).raw();
    };

    boost::numeric::odeint::integrate_times(
            Stepper{}, rhs, initialState,
            tEval.begin(), tEval.end(), dt, observer);
    return result;
}

}  // namespace country
}  // namespace epidemics
