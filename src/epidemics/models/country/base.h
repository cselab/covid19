#pragma once

#include <epidemics/integrator.h>

#include <boost/array.hpp>
#include <vector>

namespace epidemics {
namespace country {

/** Design parameters, shared between all models.
 *
 * As opposed to model parameters, design parameters are not fitted.
 */
struct DesignParameters {
    int N;  // Country population.
};


/** A base of State classes.
 *
 * The reason it exists is because setting up a custom type to work with
 * Boost's ODE integrator was not simple. As a workaround, our State class
 * (which has a nice .S(), .I(), .R()... interface) is a wrapper around a
 * boost::array (RawState).
 */
template <typename T, size_t N>
struct StateBase {
    using RawState = boost::array<T, N>;

    explicit StateBase() { }
    explicit StateBase(RawState state) : v_{std::move(state)} { }

    const RawState &raw() const noexcept { return v_; }
    RawState& raw() noexcept { return v_; }

    static constexpr size_t size() noexcept { return N; }

protected:
    RawState v_;
};


/** CRTP base class for solvers.
 *
 * Solvers have to only define a `rhs` function, the integrator is handled by
 * the base class in `base.hh`.
 */
template <typename Derived,
          template <typename> class State,
          template <typename> class Parameters>
class SolverBase {
public:
    SolverBase(DesignParameters dp) : dp_{std::move(dp)} { }

    const DesignParameters &designParameters() const noexcept { return dp_; }

    template <typename T>
    std::vector<State<T>> solve(
            const Parameters<T> &parameters,
            State<T> y0,
            const std::vector<double> &tEval,
            IntegratorSettings settings) const
    {
        return integrate(
                [this, parameters](double t, const State<T> &x, State<T> &dxdt) {
                    return derived()->rhs(t, parameters, x, dxdt);
                },
                std::move(y0), tEval, std::move(settings));
    }

protected:
    Derived *derived() noexcept {
        return static_cast<Derived *>(this);
    }
    const Derived *derived() const noexcept {
        return static_cast<const Derived *>(this);
    }

    DesignParameters dp_;
};

}  // namespace country
}  // namespace epidemics
