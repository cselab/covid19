#pragma once

#include <epidemics/integrator.h>

#include <boost/array.hpp>
#include <vector>

namespace epidemics {
namespace country {

/** Model data, shared between all models.
 *
 * Model data is not the same as model parameters!
 * These values are not fitted.
 */
struct ModelData {
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

    const RawState &raw() const { return v_; }
    RawState& raw() { return v_; }

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
          template <typename> typename State,
          template <typename> typename Parameters>
class SolverBase {
public:
    SolverBase(ModelData data) : data_{std::move(data)} { }

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

    ModelData data_;
};

}  // namespace country
}  // namespace epidemics
