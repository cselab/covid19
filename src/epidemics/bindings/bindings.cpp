#include "bindings.h"

#include <epidemics/models/country/base.h>
#include <epidemics/models/cantons/data.h>
#include <epidemics/integrator.h>

namespace py = pybind11;
using namespace py::literals;

namespace epidemics {

IntegratorSettings integratorSettingsFromKwargs(py::kwargs kwargs) {
    auto pop = kwargs.attr("pop");
    IntegratorSettings out;
    out.dt = pop("dt", out.dt).cast<double>();
    if (!kwargs.empty())
        throw py::key_error(kwargs.begin()->first.cast<std::string>());
    return out;
}


namespace country {

void exportCountryModels(py::module &top, py::module &m);

static void exportDesignParameters(py::module &m) {
    py::class_<DesignParameters>(m, "DesignParameters")
        .def(py::init<int>(), "N"_a)
        .def_readwrite("N", &DesignParameters::N, "Country population.");
}

}  // namespace country


namespace cantons {

void exportCantonsModels(py::module &top, py::module &m);

static void exportDesignParameters(py::module &m) {
    py::class_<DesignParameters>(m, "DesignParameters")
        .def(py::init<std::vector<std::string>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>>(),
             "region_keys"_a, "Ni"_a, "Mij"_a, "Cij"_a,
             "ext_com_iu"_a, "Ui"_a)
        .def_readonly("Mij", &DesignParameters::Mij)
        .def_readonly("num_regions", &DesignParameters::numRegions);
}

}  // namespace cantons
}  // namespace epidemics

PYBIND11_MODULE(libepidemics, m)
{
    using namespace epidemics;
    py::class_<IntegratorSettings>(m, "IntegratorSettings")
        .def(py::init<double>(), "dt"_a)
        .def_readwrite("dt", &IntegratorSettings::dt);

    auto country = m.def_submodule("country");
    epidemics::country::exportCountryModels(m, country);
    epidemics::country::exportDesignParameters(country);

    auto cantons = m.def_submodule("cantons");
    epidemics::cantons::exportCantonsModels(m, cantons);
    epidemics::cantons::exportDesignParameters(cantons);
}
