#pragma once

#include <epidemics/utils/autodiff.h>

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
template <typename T, int N>
struct StateBase {
    using RawState = boost::array<T, N>;

    explicit StateBase(RawState state) :
        v_{std::move(state)}
    { }

    const RawState &raw() const & { return v_; }
    RawState raw() && { return std::move(v_); }

protected:
    RawState v_;
};


/** CRTP base class for solvers.
 *
 * Solvers have to only define a `rhs` function, the integrator is handled by
 * the base class in `base.hh`.
 */
template <typename Derived, typename State, typename Parameters>
class SolverBase {
public:
    SolverBase(ModelData data) : data_{std::move(data)} { }

    std::vector<State> solve(const Parameters &parameters,
                             typename State::RawState initialState,
                             const std::vector<double> &tEval) const;
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
