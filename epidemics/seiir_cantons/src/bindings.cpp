#include "model.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(libseiir, m)
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

    py::class_<MultiSEIIR>(m, "MultiSEIIR")
        .def(py::init<std::vector<double>>(),
             "Mij"_a)
        .def("solve", &MultiSEIIR::solve,
             "parameters"_a, "initialState"_a, "days"_a);
}
