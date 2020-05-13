/// Export via pybind11 everything in `libepidemics.country.{{NAME}}`.

#include "bindings.h"
#include "country.h"
#include <epidemics/models/country/base.hh>
#include <epidemics/models/country/{{NAME}}.h>

using namespace py::literals;

namespace epidemics {
namespace country {
namespace {{NAME}} {
    template <typename T>
    static auto exportState(py::module &m, const char *name) {
        return exportGenericState<State<T>>(m, name)
        {%- for field in STATE %}
            .def("{{field}}", py::overload_cast<>(&State<T>::{{field}}, py::const_), "Get {{field}}.")
        {%- endfor -%}
        ;
    }

    template <typename T>
    static auto exportParameters(py::module &m, const char *name) {
        return py::class_<Parameters<T>>(m, name)
            .def(py::init<{{ (['T'] * PARAMS|length) | join(', ') }}>(),
            {%- for param in PARAMS %}
                 "{{param}}"_a{% if not loop.last %},{% endif %}
            {%- endfor -%}
                )
        {%- for param in PARAMS %}
            .def_readwrite("{{param}}", &Parameters<T>::{{param}})
        {%- endfor -%}
        ;
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
            }), "params"_a, "Create an AD parameters struct from non-ad parameters.");

        exportState<double>(m, "State");
        exportState<AD>(m, "StateAD")
            .def(py::init(&convertScalarStateToAD<State, AD>), "state"_a,
                 "Convert scalar state to AD state");

        py::implicitly_convertible<Parameters<double>, Parameters<AD>>();
        py::implicitly_convertible<State<double>, State<AD>>();

        m.attr("ElementAD") = exportAutoDiff<AD>(top);
        exportSolver<Solver, State, Parameters>(m);
    }
}  // namespace {{NAME}}
}  // namespace country
}  // namespace epidemics

