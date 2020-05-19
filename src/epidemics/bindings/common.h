/// Code common to country and cantons models.

#pragma once

#include "bindings.h"
#include <epidemics/integrator.hh>

namespace epidemics {

/// Shorthand for the autodiff type used for the given model.
template <template <typename> class Parameters>
using ADType = AutoDiff<double, Parameters<double>::numParameters>;

IntegratorSettings integratorSettingsFromKwargs(py::kwargs kwargs);

template <typename Solver, typename State, typename Parameters, typename PySolver>
void exportSolverSolve(
        PySolver &pySolver,
        const char *name)
{
    using namespace py::literals;
    pySolver.def(
            name,
            [](const Solver &solver,
               const Parameters &params,
               State state,
               const std::vector<double> &tEval,
               py::kwargs kwargs)
            {
                SignalRAII breakRAII;
                return solver.solve(params, std::move(state), tEval,
                                    integratorSettingsFromKwargs(kwargs));
            }, "parameters"_a, "initial_state"_a, "t_eval"_a);
}

/// Export a model Solver. Returns the Solver class handler.
template <typename Solver,
          typename ModelData,
          template <typename> class State,
          template <typename> class Parameters>
auto exportSolver(py::module &m) {
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
        }, "Convert to a Python list of elements.");
}

}  // namespace epidemics
