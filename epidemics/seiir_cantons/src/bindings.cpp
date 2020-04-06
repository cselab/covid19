#include "model.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(libseiir, m)
{
    using namespace pybind11::literals;
    m.def("sir_infected_so_far", &sir_infected_so_far,
          "S0"_a, "I0"_a, "R0"_a, "days"_a, "N"_a, "beta"_a, "gamma"_a);
}
