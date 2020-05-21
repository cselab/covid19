#pragma once

#include "bindings.h"
#include <epidemics/utils/autodiff.h>

#include <sstream>

namespace epidemics {

/// Export AutoDiff<T, N> to Python.
template <typename AD>
py::handle exportAutoDiff(py::module &m) {
    static py::handle cls;
    if (cls)
        return cls;

    using namespace py::literals;
    using T = typename AD::ValueType;
    constexpr int N = AD::kNumVariables;

    // FIXME: Name collision will occur if multiple different Ts are used.
    std::string name = "AutoDiff_T_" + std::to_string(N);
    cls = py::class_<AD>(m, name.c_str())
        .def(py::init<T>(), "val"_a, "AutoDiff with all derivatives set to 0.")
        .def(py::init<T, std::array<T, N>>(), "val"_a, "d"_a)
        .def("__repr__", [name](const AD &ad) {
            std::stringstream out;
            out << name << '(' << ad.val();
            for (int i = 0; i < N; ++i)
                out << (i == 0 ? "; " : ", ") << ad.d(i);
            out << ')';
            return std::move(out).str();
        })
        .def("val", py::overload_cast<>(&AD::val, py::const_), "Get value.")
        .def("d", py::overload_cast<int>(&AD::d, py::const_), "i"_a,
             "Get the derivative with respect to the ith variable.")
        .def("d", [](const AD &ad) {
                return std::vector<T>(ad.dBegin(), ad.dEnd());
        }, "Return a list of all derivatives.");

    py::implicitly_convertible<long long, AD>();
    py::implicitly_convertible<T, AD>();
    return cls;
}

}  // namespace epidemics
