#include <epidemics/models/country/base.hh>
#include <epidemics/models/country/sir.h>

#include "bindings.h"
#include "autodiff.h"

using namespace py::literals;

namespace epidemics {
namespace country {

/// Export everything in `libepidemics.country.sir`.
static void exportSIR(py::module &top, py::module &m) {
    using namespace sir;
    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double>(), "beta"_a, "gamma"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("gamma", &Parameters::gamma);

    m.attr("Element") = exportAutoDiff<sir::State::RawState::value_type>(top);
    py::class_<State>(m, "State")
        .def(py::init<State::RawState>())
        .def("tolist", [](const State &state) {
            return state.raw();
        }, "Convert to a Python list of floats.")
        .def("S", py::overload_cast<>(&State::S, py::const_), "Get S.")
        .def("I", py::overload_cast<>(&State::I, py::const_), "Get I.")
        .def("R", py::overload_cast<>(&State::R, py::const_), "Get R.");

    py::class_<Solver>(m, "Solver")
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


/// Export all country models to Python.
void exportCountryModels(py::module &top, py::module &m)
{
    py::class_<ModelData>(m, "ModelData")
        .def(py::init<int>(), "N"_a)
        .def_readwrite("N", &ModelData::N, "Country population.");

    auto sir = m.def_submodule("sir");
    exportSIR(top, sir);
}


}  // namespace country
}  // namespace epidemics
