#include "bindings.h"

namespace py = pybind11;

void exportCantonModels(py::module &m);

PYBIND11_MODULE(libepidemics, m)
{
    auto cantons = m.def_submodule("cantons");
    exportCantonModels(cantons);
}
