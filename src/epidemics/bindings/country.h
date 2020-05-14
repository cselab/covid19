#include "autodiff.h"
#include "common.h"

#include <epidemics/models/country/base.h>

namespace epidemics {
namespace country {

template <template <typename> typename Parameters>
using ADType = AutoDiff<double, Parameters<double>::numParameters>;

/// Export a model Solver. Returns the Solver class handler.
template <typename Solver,
          template <typename> typename State,
          template <typename> typename Parameters>
static auto exportSolver(py::module &m) {
    using namespace py::literals;
    using AD = ADType<Parameters>;

    auto solver = py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData>(), "model_data"_a);
    exportSolverSolve<Solver, State<double>, Parameters<double>>(solver, "solve");
    exportSolverSolve<Solver, State<AD>, Parameters<AD>>(solver, "solve_ad");
    return solver;
}


/// Export a model State. Returns the State class handler.
template <typename State>
static auto exportGenericState(py::module &m, const char *name) {
    return py::class_<State>(m, name)
        .def(py::init<typename State::RawState>())
        .def("tolist", [](const State &state) {
            return state.raw();
        }, "Convert to a Python list of floats.");
}

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
