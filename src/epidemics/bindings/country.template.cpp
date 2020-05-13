#include "bindings.h"
#include <epidemics/models/country/base.h>

using namespace py::literals;

namespace epidemics {
namespace country {

{% for model in MODELS %}
namespace {{model}} { void exportAll(py::module &top, py::module &m); }
{% endfor %}

/// Export all country models to Python.
void exportCountryModels(py::module &top, py::module &m)
{
    py::class_<ModelData>(m, "ModelData")
        .def(py::init<int>(), "N"_a)
        .def_readwrite("N", &ModelData::N, "Country population.");

    {% for model in MODELS %}
    auto {{model}} = m.def_submodule("{{model}}");
    {{model}}::exportAll(top, {{model}});
    {% endfor %}
}

}  // namespace country
}  // namespace epidemics

