#include <epidemics/models/country/base.hh>
#include <epidemics/models/country/seiir.h>
#include <epidemics/models/country/sir.h>

#include "bindings.h"
#include "autodiff.h"

using namespace py::literals;

namespace epidemics {
namespace country {


/// Export a model Solver. Returns the Solver class handler.
template <typename Solver, typename State, typename Parameters>
static auto exportSolver(py::module &m) {
    return py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData>(), "model_data"_a)
        .def("solve", [](const Solver &solver,
                         const Parameters &params,
                         State state,
                         const std::vector<double> &tEval)
            {
                SignalRAII break_raii;
                return solver.solve(params, std::move(state).raw(), tEval);
            },
            "parameters"_a, "initial_state"_a, "t_eval"_a);
}


/// Export a model State. Returns the State class handler.
template <typename State>
static auto exportState(py::module &m) {
    return py::class_<State>(m, "State")
        .def(py::init<typename State::RawState>())
        .def("tolist", [](const State &state) {
            return state.raw();
        }, "Convert to a Python list of floats.");
}

/// Export everything in `libepidemics.country.sir`.
static void exportSIR(py::module &top, py::module &m) {
    using namespace sir;
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double>(), "beta"_a, "gamma"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("gamma", &Parameters::gamma);

    m.attr("Element") = exportAutoDiff<State::RawState::value_type>(top);
    exportState<State>(m)
        .def("S", py::overload_cast<>(&State::S, py::const_), "Get S.")
        .def("I", py::overload_cast<>(&State::I, py::const_), "Get I.")
        .def("R", py::overload_cast<>(&State::R, py::const_), "Get R.");

    exportSolver<Solver, State, Parameters>(m);
}

/// Export everything in `libepidemics.country.seiir`.
static void exportSEIIR(py::module &top, py::module &m) {
    using namespace seiir;
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double, double, double, double>(),
             "beta"_a, "mu"_a, "alpha"_a, "Z"_a, "D"_a)
        .def_readwrite("beta",  &Parameters::beta)
        .def_readwrite("mu",    &Parameters::mu)
        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("Z",     &Parameters::Z)
        .def_readwrite("D",     &Parameters::D);

    m.attr("Element") = exportAutoDiff<State::RawState::value_type>(top);
    exportState<State>(m)
        .def("S",  py::overload_cast<>(&State::S,  py::const_), "Get S.")
        .def("E",  py::overload_cast<>(&State::E,  py::const_), "Get E.")
        .def("Ir", py::overload_cast<>(&State::Ir, py::const_), "Get Ir.")
        .def("Iu", py::overload_cast<>(&State::Iu, py::const_), "Get Iu.")
        .def("R",  py::overload_cast<>(&State::R,  py::const_), "Get R.");

    exportSolver<Solver, State, Parameters>(m);
}

/// Export all country models to Python.
void exportCountryModels(py::module &top, py::module &m)
{
    py::class_<ModelData>(m, "ModelData")
        .def(py::init<int>(), "N"_a)
        .def_readwrite("N", &ModelData::N, "Country population.");

    auto sir   = m.def_submodule("sir");
    auto seiir = m.def_submodule("seiir");
    exportSIR  (top, sir);
    exportSEIIR(top, seiir);
}


}  // namespace country
}  // namespace epidemics
