/// Code common to country and cantons models.

#pragma once

#include "bindings.h"
#include <epidemics/integrator.hh>

namespace epidemics {

/// Shorthand for the autodiff type used for the given model.
template <template <typename> class Parameters>
using StaticADType = AutoDiff<double, Parameters<double>::numParameters>;

IntegratorSettings integratorSettingsFromKwargs(py::kwargs kwargs);

/// Attributes common to both non-AD and AD variants.
template <typename Solver, typename State, typename Parameters, typename PySolver>
void exportSolverCommon(
        py::module &m,
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
            }, "params"_a, "y0"_a, "t_eval"_a);
    pySolver.attr("model") = m;
}

/// Export a model Solver. Returns the Solver class handler.
template <typename Solver,
          typename DesignParameters,
          template <typename> class State,
          template <typename> class Parameters>
auto exportSolver(py::module &m) {
    using namespace py::literals;
    using StaticAD = StaticADType<Parameters>;
    using DynamicAD = DynamicAutoDiff<double>;

    auto solver = py::class_<Solver>(m, "Solver")
        .def(py::init<DesignParameters>(), "dp"_a)
        .def_property_readonly("dp", &Solver::designParameters)
        .def("solve_ad", [](const Solver &, py::args, py::kwargs) {
            throw std::runtime_error("'solve_ad' has been renamed to 'solve_params_ad'");
        });
    exportSolverCommon<Solver, State<double>, Parameters<double>>(m, solver, "solve");
    exportSolverCommon<Solver, State<StaticAD>, Parameters<StaticAD>>(m, solver, "solve_params_ad");

    // You shouldn't use _solve_custom_ad directly, use {country,cantons}_custom_derivatives instead.
    exportSolverCommon<Solver, State<DynamicAD>, Parameters<DynamicAD>>(m, solver, "_solve_custom_ad");
    return solver;
}

/// Export a model State. Returns the State class handler.
template <typename State>
static auto exportGenericState(py::module &m, const char *name) {
    return py::class_<State>(m, name)
        .def(py::init<typename State::RawState>())
        .def("tolist", [](const State &state) {
            return state.raw();
        }, "Convert to a Python list of elements.")
        .def("__call__", [](const State &state, size_t index) {
            if (index < state.raw().size())
                return state.raw()[index];
            throw std::out_of_range(std::to_string(index));
        }, "Get state variables by index. We don't use __getitem__, because "
           "it would break the implicit cast from AD states to non-AD states.")
        .def("__len__", &State::size, "Get the total number of state variables.");
}

}  // namespace epidemics
