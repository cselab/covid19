#pragma once

#include "data.h"

#include <cassert>

namespace epidemics {
namespace cantons {

/*
 * Using custom types with boost::odeint is not that simple. Instead, we use a
 * single std::vector<> to store the whole state of the simulation, and the
 * wrapper below to access the data in a human-friendly way.
 */
template <typename T, size_t VarsPerRegion>
struct StateBase {
    using RawState = std::vector<T>;
    static constexpr size_t kVarsPerRegion = VarsPerRegion;

    /// Create a state with all values set to 0.
    explicit StateBase(size_t numRegions) :
        numRegions_{numRegions},
        v_(kVarsPerRegion * numRegions, 0.0)
    { }

    /// Create a state from a raw state.
    explicit StateBase(RawState state) :
        numRegions_{state.size() / kVarsPerRegion},
        v_{std::move(state)}
    {
        assert(v_.size() % kVarsPerRegion == 0);
    }

    size_t numRegions() const noexcept { return numRegions_; }
    RawState &raw() & { return v_; }
    const RawState &raw() const & { return v_; }
    RawState raw() && { return std::move(v_); }

protected:
    size_t numRegions_;
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
    SolverBase(ModelData modelData) :
        modelData_{std::move(modelData)}
    { }

    double M(int from, int to) const {
        return modelData_.Mij[from * modelData_.numRegions + to];
    }
    const std::vector<size_t>& nonzero_Mij(size_t from_or_in) const {
        return modelData_.nonzero_Mij[from_or_in];
    }
    double C(int from, int to) const {
        return modelData_.Cij[from * modelData_.numRegions + to];
    }
    double C_plus_Ct(int from, int to) const {
        return modelData_.C_plus_Ct[from * modelData_.numRegions + to];
    }

    template <typename T>
    std::vector<State<T>> solve(const Parameters<T> &parameters,
                                typename State<T>::RawState initialState,
                                const std::vector<double> &tEval) const;

protected:
    Derived *derived() noexcept {
        return static_cast<Derived *>(this);
    }
    const Derived *derived() const noexcept {
        return static_cast<const Derived *>(this);
    }

    ModelData modelData_;
};

}  // namespace cantons
}  // namespace epidemics
