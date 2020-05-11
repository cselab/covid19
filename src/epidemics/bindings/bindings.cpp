#include "bindings.h"

namespace epidemics {

namespace py = pybind11;

namespace country {
    void exportCountryModels(py::module &top, py::module &m);
};
namespace cantons {
    void exportCantonModels(py::module &top, py::module &m);
};

}  // namespace epidemics

PYBIND11_MODULE(libepidemics, m)
{
    auto cantons = m.def_submodule("cantons");
    epidemics::cantons::exportCantonModels(m, cantons);

    auto country = m.def_submodule("country");
    epidemics::country::exportCountryModels(m, country);
}
