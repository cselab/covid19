#include "autodiff.h"
#include "common.h"

#include <epidemics/models/country/base.h>

namespace epidemics {
namespace country {

template <template <typename> typename State, typename AD>
State<AD> convertScalarStateToAD(const State<double> &state) {
    State<AD> out;
    static_assert(out.size() == state.size(), "Sanity check.");
    for (size_t i = 0; i < state.size(); ++i)
        out.raw()[i] = state.raw()[i];
    return out;
}

}  // namespace country
}  // namespace epidemics
