#pragma once

#include <epidemics/utils/autodiff.h>
#include <epidemics/utils/signal.h>
#include <boost/numeric/odeint.hpp>

#include <vector>

namespace epidemics {

struct IntegratorSettings {
    double dt{0.1};
};

template <typename RHS, typename State>
std::vector<State> integrate(
        RHS rhs,
        State y0,
        const std::vector<double> &tEval,
        IntegratorSettings settings);
}  // namespace epidemics


namespace boost {
namespace numeric {
namespace odeint {

/// Partial template specialization of boost's internal vector resize functor
/// for DynamicAD types, which need to know their size at construct time.
template <class T>
#if BOOST_VERSION >= 105600
struct resize_impl_sfinae
#else
struct resize_impl
#endif
    <std::vector<epidemics::DynamicAutoDiff<T>>,
     std::vector<epidemics::DynamicAutoDiff<T>>>
{
    using AD = epidemics::DynamicAutoDiff<T>;
    using State = std::vector<AD>;
    static void resize(State &x1, const State &x2)
    {
        x1.resize(x2.size(), AD((typename AD::size_tag){}, x2[0].N()));
    }
};
}  // namespace odeint
}  // namespace numeric
}  // namespace boost
