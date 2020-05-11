#include "bindings.h"

namespace py = pybind11;

void exportCantonModels(py::module &m);

PYBIND11_MODULE(libsolver, m)
{
    // TODO: Restructure namespaces.
    exportCantonModels(m);
}
