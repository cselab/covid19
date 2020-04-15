#include "model.h"
#include "solver.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

auto makeValuesGetter(int valueIndex) {
    assert(0 <= valueIndex && valueIndex < kVarsPerRegion);
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

PYBIND11_MODULE(libsolver, m)
{
    using namespace pybind11::literals;

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<double, double, double, double, double, double>(),
             "beta"_a, "mu"_a, "alpha"_a, "Z"_a, "D"_a, "theta"_a)
        .def_readwrite("beta", &Parameters::beta)
        .def_readwrite("mu", &Parameters::mu)
        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("Z", &Parameters::Z)
        .def_readwrite("D", &Parameters::D)
        .def_readwrite("theta", &Parameters::theta);

    py::class_<ModelData>(m, "ModelData")
        .def(py::init<size_t, std::map<std::string, int>, std::vector<int>,
                      std::vector<double>, std::vector<double>>());
    m.def("readModelData", &readModelData, "filename");

    py::class_<State>(m, "State")
        .def(py::init<std::vector<double>>())
        .def("tolist", [](const State &state) -> std::vector<double> {
            return state.raw();
        }, "Convert to a Python list of floats.")
        .def("S", makeValuesGetter(0), "Get a list of S for each region.")
        .def("E", makeValuesGetter(1), "Get a list of E for each region.")
        .def("Ir", makeValuesGetter(2), "Get a list of Ir for each region.")
        .def("Iu", makeValuesGetter(3), "Get a list of Iu for each region.")
        .def("N", makeValuesGetter(4), "Get a list of N for each region.")
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
