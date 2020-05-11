#include "bindings.h"

namespace epidemics {

namespace py = pybind11;

namespace cantons { void exportCantonModels(py::module &m); };

}  // namespace epidemics

PYBIND11_MODULE(libepidemics, m)
{
    auto cantons = m.def_submodule("cantons");
    epidemics::cantons::exportCantonModels(cantons);
}
