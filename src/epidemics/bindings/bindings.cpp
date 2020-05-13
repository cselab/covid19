#include "bindings.h"

#include <epidemics/models/country/base.h>
#include <epidemics/models/cantons/data.h>

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

void exportCantonsModels(py::module &top, py::module &m);

static void exportModelData(py::module &m) {
    py::class_<ModelData>(m, "ModelData")
        .def(py::init<std::vector<std::string>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>>(),
             "region_keys"_a, "Ni"_a, "Mij"_a, "Cij"_a,
             "ext_com_iu"_a, "Ui"_a)
        .def_readonly("Mij", &ModelData::Mij)
        .def_readonly("num_regions", &ModelData::numRegions);
}

}  // namespace cantons
}  // namespace epidemics

PYBIND11_MODULE(libepidemics, m)
{
    auto country = m.def_submodule("country");
    epidemics::country::exportCountryModels(m, country);
    epidemics::country::exportModelData(country);

    auto cantons = m.def_submodule("cantons");
    epidemics::cantons::exportCantonsModels(m, cantons);
    epidemics::cantons::exportModelData(cantons);
}
