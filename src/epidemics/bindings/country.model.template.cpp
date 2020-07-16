/// Export via pybind11 everything in `libepidemics.country.{{NAME}}`.

#include "bindings.h"
#include "country.h"
#include <epidemics/models/country/base.h>
#include <epidemics/models/country/{{NAME}}.h>

using namespace py::literals;

namespace epidemics {
namespace country {
namespace {{NAME}} {

{{EXPORT_PARAMETERS}}

template <typename T>
static auto exportState(py::module &m, const char *name) {
    return exportGenericState<State<T>>(m, name)
    {%- for field in STATE %}
        .def("{{field}}", py::overload_cast<>(&State<T>::{{field}}, py::const_), "Get {{field}}.")
    {%- endfor -%}
    ;
}

void exportAll(py::module &top, py::module &m) {
    using StaticAD = StaticADType<Parameters>;
    using DynamicAD = DynamicAutoDiff<double>;
    exportParameters<double>(m, "Parameters");
    exportParameters<StaticAD>(m, "ParametersStaticAD")
        .def(py::init([](const Parameters<double> &params) {
            return Parameters<StaticAD>{
                    {% for field in PARAMS %}
                        {%- set outer_loop = loop -%}
                        make_ad(params.{{field}}
                        {%- for inner in PARAMS -%}
                            {% if outer_loop.index == loop.index %}, 1{% else %}, 0{% endif -%}
                        {%- endfor -%}),
                    {% endfor %}
            };
        }), "params"_a, "Create a StaticAD parameters struct from non-ad parameters.");
    exportParameters<DynamicAD>(m, "ParametersDynamicAD")
        .def(py::init([](const Parameters<double> &) -> Parameters<DynamicAD> {
            // TODO: Is this even necessary?
            throw std::invalid_argument("Parameters -> ParametersDynamicAD not yet implemented.");
        }), "params"_a, "Create a DynamicAD parameters struct from non-ad parameters.");

    exportState<double>(m, "State");
    exportState<StaticAD>(m, "StateStaticAD")
        .def(py::init(&convertScalarStateToAD<State, StaticAD>), "state"_a,
             "Convert scalar state to StaticAD state");
    exportState<DynamicAD>(m, "StateDynamicAD");

    py::implicitly_convertible<Parameters<double>, Parameters<StaticAD>>();
    py::implicitly_convertible<State<double>, State<StaticAD>>();

    m.attr("StaticAD") = exportStaticAutoDiff<StaticAD>(top, "StaticAD_double_");
    m.attr("DynamicAD") = exportDynamicAutoDiff<DynamicAD>(top, "DynamicAD_double_");
    exportSolver<Solver, DesignParameters, State, Parameters>(m)
        .def("state_size", [](const Solver &) noexcept {
            return State<double>::size();
        }, "Return the number of state variables.");
}

}  // namespace {{NAME}}
}  // namespace country
}  // namespace epidemics
