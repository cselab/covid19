#include "bindings.h"
#include "autodiff.h"

#include <epidemics/models/cantons/data.h>

namespace epidemics {
namespace cantons {

template <template <typename> typename Parameters>
using ADType = AutoDiff<double, Parameters<double>::numParameters>;

/// Export a model Solver. Returns the Solver class handler.
template <typename Solver,
          template <typename> typename State,
          template <typename> typename Parameters>
static auto exportSolver(py::module &m) {
    using namespace py::literals;
    using AD = ADType<Parameters>;

    return py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData>(), "model_data"_a)
        .def("solve", [](const Solver &solver,
                         const Parameters<double> &params,
                         State<double> state,
                         const std::vector<double> &tEval)
            {
                SignalRAII breakRAII;
                return solver.solve(params, std::move(state).raw(), tEval);
            },
            "parameters"_a, "initial_state"_a, "t_eval"_a)
        .def("solve_ad", [](const Solver &solver,
                            const Parameters<AD> &params,
                            State<AD> state,
                            const std::vector<double> &tEval)
            {
                SignalRAII breakRAII;
                return solver.solve(params, std::move(state).raw(), tEval);
            },
            "parameters"_a, "initial_state"_a, "t_eval"_a);
}


/// Helper function factory for State getters.
template <template <typename> typename State, typename T>
static auto makeValuesGetter(size_t valueIndex) {
    assert(0 <= valueIndex && valueIndex < State<T>::kVarsPerRegion);
    /// Extract a subvector of the state corresponding to the given value.
    return [valueIndex](const State<T> &state) {
        const T *p = state.raw().data();
        return std::vector<T>(p + valueIndex * state.numRegions(),
                              p + (valueIndex + 1) * state.numRegions());
    };
}

/// Export a model State. Returns the State class handler.
template <typename State>
static auto exportGenericState(py::module &m, const char *name) {
    return py::class_<State>(m, name)
        .def(py::init<typename State::RawState>())
        .def("tolist", [](const State &state) {
            return state.raw();
        }, "Convert to a Python list of elements.");
}

template <template <typename> typename State, typename AD>
State<AD> convertScalarStateToAD(const State<double> &state) {
    State<AD> out(state.numRegions());
    for (size_t i = 0; i < state.raw().size(); ++i)
        out.raw()[i] = state.raw()[i];
    return out;
}

}  // namespace cantons
}  // namespace epidemics
