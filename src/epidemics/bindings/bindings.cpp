#include "bindings.h"

#include <epidemics/models/country/base.h>

namespace py = pybind11;
using namespace py::literals;

namespace epidemics {
namespace country {

void exportCountryModels(py::module &top, py::module &m);

static void exportModelData(py::module &m) {
    py::class_<ModelData>(m, "ModelData")
        .def(py::init<int>(), "N"_a)
        .def_readwrite("N", &ModelData::N, "Country population.");
}

}  // namespace country

namespace cantons {
    void exportCantonModels(py::module &top, py::module &m);
}  // namespace cantons

}  // namespace epidemics

PYBIND11_MODULE(libepidemics, m)
{
    auto cantons = m.def_submodule("cantons");
    epidemics::cantons::exportCantonModels(m, cantons);

    auto country = m.def_submodule("country");
    epidemics::country::exportCountryModels(m, country);
    epidemics::country::exportModelData(country);
}
