/// Export via pybind11 everything in `libepidemics.cantons.{{NAME}}`.

#include "bindings.h"
#include "cantons.h"
#include <epidemics/models/cantons/base.h>
#include <epidemics/models/cantons/{{NAME}}.h>

using namespace py::literals;

namespace epidemics {
namespace cantons {
namespace {{NAME}} {

{{EXPORT_PARAMETERS}}

template <typename T>
static auto exportState(py::module &m, const char *name) {
    return exportGenericState<State<T>>(m, name)
    {%- for field in STATE %}
        .def("{{field}}", makeValuesGetter<State, T>({{ loop.index0 }}), "Get a list of {{field}} for each canton.")
    {%- endfor -%}
    {%- for field in STATE %}
        .def("{{field}}", py::overload_cast<size_t>(&State<T>::{{field}}, py::const_), "Get {{field}}_i.")
    {%- endfor -%};
}

void exportAll(py::module &top, py::module &m) {
    using AD = ADType<Parameters>;
    exportParameters<double>(m, "Parameters");
    exportParameters<AD>(m, "ParametersAD")
        .def(py::init([](const Parameters<double> &params) {
            return Parameters<AD>{
                    {% for field in PARAMS %}
                    {%- set outer_loop = loop -%}
                    make_ad(params.{{field}}
                        {%- for inner in PARAMS -%}
                            {% if outer_loop.index == loop.index %}, 1{% else %}, 0{% endif -%}
                        {%- endfor -%}),
                    {% endfor %}
            };
        }), "params"_a, "Create an AD parameters struct from non-AD parameters.");

    exportState<double>(m, "State");
    exportState<AD>(m, "StateAD")
        .def(py::init(&convertScalarStateToAD<State, AD>), "state"_a,
             "Convert scalar state to an AD state");

    py::implicitly_convertible<Parameters<double>, Parameters<AD>>();
    py::implicitly_convertible<State<double>, State<AD>>();

    m.attr("AD") = exportAutoDiff<AD>(top);
    exportSolver<Solver, ModelData, State, Parameters>(m)
        .def("state_size", &Solver::stateSize, "Return the number of state variables.");
}

}  // namespace {{NAME}}
}  // namespace cantons
}  // namespace epidemics
