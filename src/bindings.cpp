#include "common.hh"
#include "data.h"
#include "model_seiin.h"
#include "model_seii_c.h"
#include "model_sei_c.h"
#include "model_seiin_interventions.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

template <typename State>
auto makeValuesGetter(int valueIndex) {
    assert(0 <= valueIndex && valueIndex < State::kVarsPerRegion);
    /// Extract a subvector of the state corresponding to the given value.
    return [valueIndex](const State &state) {
        const double *p = state.raw().data();
        return std::vector<double>(p + valueIndex * state.numRegions(),
                                   p + (valueIndex + 1) * state.numRegions());
    };
}

class SignalRAII {
public:
    SignalRAII() {
        check_signals_func = []() {
            // https://stackoverflow.com/questions/14707049/allowing-ctrl-c-to-interrupt-a-python-c-extension
            if (PyErr_CheckSignals() != 0)
                throw std::runtime_error("Signal received. Breaking.");
        };
    }
    ~SignalRAII() {
        check_signals_func = nullptr;
    }
};

void exportSEIIN(py::module &m) {
    using namespace pybind11::literals;
    using namespace seiin;

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double, double, double, double, double>(),
             "beta"_a, "mu"_a, "alpha"_a, "Z"_a, "D"_a, "theta"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("mu", &Parameters::mu)
        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("Z", &Parameters::Z)
        .def_readwrite("D", &Parameters::D)
        .def_readwrite("theta", &Parameters::theta);

    py::class_<State>(m, "State")
        .def(py::init<std::vector<double>>())
        .def("tolist", [](const State &state) -> std::vector<double> {
            return state.raw();
        }, "Convert to a Python list of floats.")
        .def("S", makeValuesGetter<State>(0), "Get a list of S for each region.")
        .def("E", makeValuesGetter<State>(1), "Get a list of E for each region.")
        .def("Ir", makeValuesGetter<State>(2), "Get a list of Ir for each region.")
        .def("Iu", makeValuesGetter<State>(3), "Get a list of Iu for each region.")
        .def("N", makeValuesGetter<State>(4), "Get a list of N for each region.")
        .def("S", py::overload_cast<size_t>(&State::S, py::const_), "Get S_i.")
        .def("E", py::overload_cast<size_t>(&State::E, py::const_), "Get E_i.")
        .def("Ir", py::overload_cast<size_t>(&State::Ir, py::const_), "Get Ir_i.")
        .def("Iu", py::overload_cast<size_t>(&State::Iu, py::const_), "Get Iu_i.")
        .def("N", py::overload_cast<size_t>(&State::N, py::const_), "Get N_i.");

    py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData, bool>(), "model_data"_a, "verbose"_a=false)
        .def("solve", [](const Solver &solver,
                         const Parameters &params,
                         RawState state,
                         int num_days)
            {
                SignalRAII break_raii;
                return solver.solve(params, std::move(state), num_days);
            },
            "parameters"_a, "initial_state"_a, "days"_a);
}

void exportSEII_C(py::module &m) {
    using namespace pybind11::literals;
    using namespace seii_c;

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double, double, double, double>(),
             "beta"_a, "nu"_a, "alpha"_a, "Z"_a, "D"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("nu", &Parameters::nu)
        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("Z", &Parameters::Z)
        .def_readwrite("D", &Parameters::D);

    py::class_<State>(m, "State")
        .def(py::init<std::vector<double>>())
        .def("tolist", [](const State &state) -> std::vector<double> {
            return state.raw();
        }, "Convert to a Python list of floats.")
        .def("S", makeValuesGetter<State>(0), "Get a list of S for each region.")
        .def("E", makeValuesGetter<State>(1), "Get a list of E for each region.")
        .def("Ir", makeValuesGetter<State>(2), "Get a list of Ir for each region.")
        .def("Iu", makeValuesGetter<State>(3), "Get a list of Iu for each region.")
        .def("S", py::overload_cast<size_t>(&State::S, py::const_), "Get S_i.")
        .def("E", py::overload_cast<size_t>(&State::E, py::const_), "Get E_i.")
        .def("Ir", py::overload_cast<size_t>(&State::Ir, py::const_), "Get Ir_i.")
        .def("Iu", py::overload_cast<size_t>(&State::Iu, py::const_), "Get Iu_i.");

    py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData, bool>(), "model_data"_a, "verbose"_a=false)
        .def("solve", [](const Solver &solver,
                         const Parameters &params,
                         RawState state,
                         int num_days)
            {
                SignalRAII break_raii;
                return solver.solve(params, std::move(state), num_days);
            },
            "parameters"_a, "initial_state"_a, "days"_a);
}

void exportSEI_C(py::module &m) {
    using namespace pybind11::literals;
    using namespace sei_c;

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double, double, double, double, double>(),
             "beta"_a, "nu"_a, "Z"_a, "D"_a, "tact"_a, "kbeta"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("nu", &Parameters::nu)
        .def_readwrite("Z", &Parameters::Z)
        .def_readwrite("D", &Parameters::D)
        .def_readwrite("tact", &Parameters::tact)
        .def_readwrite("kbeta", &Parameters::kbeta);

    py::class_<State>(m, "State")
        .def(py::init<std::vector<double>>())
        .def("tolist", [](const State &state) -> std::vector<double> {
            return state.raw();
        }, "Convert to a Python list of floats.")
        .def("S", makeValuesGetter<State>(0), "Get a list of S for each region.")
        .def("E", makeValuesGetter<State>(1), "Get a list of E for each region.")
        .def("I", makeValuesGetter<State>(2), "Get a list of I for each region.")
        .def("S", py::overload_cast<size_t>(&State::S, py::const_), "Get S_i.")
        .def("E", py::overload_cast<size_t>(&State::E, py::const_), "Get E_i.")
        .def("I", py::overload_cast<size_t>(&State::I, py::const_), "Get I_i.");

    py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData, bool>(), "model_data"_a, "verbose"_a=false)
        .def("solve", [](const Solver &solver,
                         const Parameters &params,
                         RawState state,
                         int num_days)
            {
                SignalRAII break_raii;
                return solver.solve(params, std::move(state), num_days);
            },
            "parameters"_a, "initial_state"_a, "days"_a);
}

void exportSEIIN_INTERVENTIONS(py::module &m) {
    using namespace pybind11::literals;
    using namespace seiin_interventions;

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double, double, double, double, double, double, double, double, double>(),
             "beta"_a, "mu"_a, "alpha"_a, "Z"_a, "D"_a, "theta"_a, "b0"_a,"b1"_a,"b2"_a,"b3"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("mu", &Parameters::mu)
        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("Z", &Parameters::Z)
        .def_readwrite("D", &Parameters::D)
        .def_readwrite("theta", &Parameters::theta)
        .def_readwrite("b0", &Parameters::b0)
        .def_readwrite("b1", &Parameters::b1)
        .def_readwrite("b2", &Parameters::b2)
        .def_readwrite("b3", &Parameters::b3);

    py::class_<State>(m, "State")
        .def(py::init<std::vector<double>>())
        .def("tolist", [](const State &state) -> std::vector<double> {
            return state.raw();
        }, "Convert to a Python list of floats.")
        .def("S", makeValuesGetter<State>(0), "Get a list of S for each region.")
        .def("E", makeValuesGetter<State>(1), "Get a list of E for each region.")
        .def("Ir", makeValuesGetter<State>(2), "Get a list of Ir for each region.")
        .def("Iu", makeValuesGetter<State>(3), "Get a list of Iu for each region.")
        .def("N", makeValuesGetter<State>(4), "Get a list of N for each region.")
        .def("S", py::overload_cast<size_t>(&State::S, py::const_), "Get S_i.")
        .def("E", py::overload_cast<size_t>(&State::E, py::const_), "Get E_i.")
        .def("Ir", py::overload_cast<size_t>(&State::Ir, py::const_), "Get Ir_i.")
        .def("Iu", py::overload_cast<size_t>(&State::Iu, py::const_), "Get Iu_i.")
        .def("N", py::overload_cast<size_t>(&State::N, py::const_), "Get N_i.");

    py::class_<Solver>(m, "Solver")
        .def(py::init<ModelData, bool>(), "model_data"_a, "verbose"_a=false)
        .def("solve", [](const Solver &solver,
                         const Parameters &params,
                         RawState state,
                         int num_days)
            {
                SignalRAII break_raii;
                return solver.solve(params, std::move(state), num_days);
            },
            "parameters"_a, "initial_state"_a, "days"_a);
}

PYBIND11_MODULE(libsolver, m)
{
    auto solvers = m.def_submodule("solvers");
    auto seiin = solvers.def_submodule("seiin");
    auto seii_c = solvers.def_submodule("seii_c");
    auto sei_c = solvers.def_submodule("sei_c");
    auto seiin_interventions = solvers.def_submodule("seiin_interventions");
    exportSEIIN(seiin);
    exportSEII_C(seii_c);
    exportSEI_C(sei_c);
    exportSEIIN_INTERVENTIONS(seiin_interventions);

    py::class_<ModelData>(m, "ModelData")
        .def(py::init<std::vector<std::string>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>>());
    m.def("readModelData", &readModelData, "filename");
}
