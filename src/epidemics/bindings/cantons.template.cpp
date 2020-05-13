#include "bindings.h"

using namespace py::literals;

namespace epidemics {
namespace cantons {

{% for model in MODELS -%}
namespace {{model}} { void exportAll(py::module &top, py::module &m); }
{% endfor %}

/// Export all cantons models to Python.
void exportCantonsModels(py::module &top, py::module &m)
{
    {% for model in MODELS %}
    auto {{model}} = m.def_submodule("{{model}}");
    {{model}}::exportAll(top, {{model}});
    {% endfor %}
}

}  // namespace cantons
}  // namespace epidemics

