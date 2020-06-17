#!/usr/bin/env python3

import jinja2
import os
PATH = os.path.normpath(os.path.dirname(__file__))

NO_EDIT_NOTE = f"""
// THIS FILE WAS CREATED FROM {os.path.basename(__file__)}
// DO NOT EDIT!

""".lstrip()

EXPORT_PARAMETERS_TEMPLATE = jinja2.Template('''
template <typename T>
static auto exportParameters(py::module &m, const char *name) {
    auto pyParams = py::class_<Parameters<T>>(m, name)
        .def(py::init<{{ (['T'] * PARAMS|length) | join(', ') }}>(),
        {%- for param in PARAMS %}
             "{{param}}"_a{% if not loop.last %},{% endif %}
        {%- endfor -%}
            )
    {%- for param in PARAMS %}
        .def_readwrite("{{param}}", &Parameters<T>::{{param}})
    {%- endfor %}
        .def("__getitem__",
            [](const Parameters<T> &params, size_t index) {
                switch (index) {
            {%- for param in PARAMS %}
                    case {{ loop.index0 }}: return params.{{param}};
            {%- endfor %}
                    default: throw py::index_error(std::to_string(index));
                }
            })
        .def("__len__", [](const Parameters<T> &) { return {{PARAMS|length}}; });
    pyParams.attr("NUM_PARAMETERS") = {{PARAMS|length}};
    return pyParams;
}
''')

with open(os.path.join(PATH, 'country.model.template.cpp')) as f:
    COUNTRY_TEMPLATE = jinja2.Template(f.read())

with open(os.path.join(PATH, 'cantons.model.template.cpp')) as f:
    CANTONS_TEMPLATE = jinja2.Template(f.read())


def save_if_modified(filename, content):
    """Save `content` to the given file, if the current content
    is different or if the file does not exist."""
    try:
        with open(filename) as f:
            current = f.read()
        if current == content:
            return
    except FileNotFoundError:
        pass

    with open(filename, 'w') as f:
        f.write(content)


def _generate_model(template, name, state, params):
    """Helper class for model generation functions."""
    state = state.split()
    params = params.split()
    kwargs = {'NAME': name, 'STATE': state, 'PARAMS': params}
    params_code = EXPORT_PARAMETERS_TEMPLATE.render(**kwargs)
    code = NO_EDIT_NOTE + template.render(
            EXPORT_PARAMETERS=params_code, **kwargs)

    # Disable because it breaks CMake's dependency checking.
    # # Make differentiation between generated and template code easier.
    # code = '/**/' + code.replace('\n', '\n/**/')
    return code


def generate_country_model(name, state, params):
    """Generate bindings for one country model."""
    code = _generate_model(COUNTRY_TEMPLATE, name, state, params)
    save_if_modified(os.path.join(PATH, f'_country.{name}.generated.cpp'), code)


def generate_canton_model(name, state, params):
    """Generate bindings for one cantons model."""
    code = _generate_model(CANTONS_TEMPLATE, name, state, params)
    save_if_modified(os.path.join(PATH, f'_cantons.{name}.generated.cpp'), code)


def generate_country(*models):
    """Generate _country.generated.cpp, responsible for libepidemics.country submodule."""
    with open(os.path.join(PATH, 'country.template.cpp')) as f:
        template = jinja2.Template(f.read())
    code = NO_EDIT_NOTE + template.render(MODELS=models)
    save_if_modified(os.path.join(PATH, '_country.generated.cpp'), code)


def generate_canton(*models):
    """Generate _canton.generated.cpp, responsible for libepidemics.cantons submodule."""
    with open(os.path.join(PATH, 'cantons.template.cpp')) as f:
        template = jinja2.Template(f.read())
    code = NO_EDIT_NOTE + template.render(MODELS=models)
    save_if_modified(os.path.join(PATH, '_cantons.generated.cpp'), code)

def main():
    # add model here
    generate_country_model('sir', 'S I R', 'beta gamma')
    generate_country_model('sir_int', 'S I R', 'beta gamma tact dtact kbeta')
    generate_country_model('sir_int_r0', 'S I R', 'r0 gamma tact dtact kbeta')
    generate_country_model('seir', 'S E I R', 'beta gamma a')
    generate_country_model('seir_int', 'S E I R', 'beta gamma a tact dtact kbeta')
    generate_country_model('seiir', 'S E Ir Iu R', 'beta mu alpha Z D')
    generate_country_model('seiir_int', 'S E Ir Iu R', 'beta mu alpha Z D tact dtact kbeta')
    # or here
    generate_canton_model('sei_c', 'S E I', 'beta nu Z D tact kbeta')
    generate_canton_model('seii_c', 'S E Ir Iu', 'beta nu alpha Z D')
    generate_canton_model('seiin', 'S E Ir Iu N', 'beta mu alpha Z D theta')
    generate_canton_model('seiin_interventions', 'S E Ir Iu N', 'beta mu alpha Z D theta b1 b2 b3 d1 d2 d3')
 
    # reparam
    generate_country_model('sir_reparam', 'S I R', 'R0 D')
    generate_country_model('sir_int_reparam', 'S I R', 'R0 D tact dtact kbeta')
    generate_country_model('sir_intexp_reparam', 'S I R', 'R0 D tact k')
    generate_country_model('sir_intrem_reparam', 'S I R', 'R0 D tact l')
    generate_country_model('seir_reparam', 'S E I R', 'R0 D Z')
    generate_country_model('seir_int_reparam', 'S E I R', 'R0 D Z tact dtact kbeta')
    generate_country_model('seir_intexp_reparam', 'S E I R', 'R0 D Z tact k')
    generate_country_model('seir_intrem_reparam', 'S E I R', 'R0 D Z tact l')
    generate_country_model('seiir_reparam', 'S E Ir Iu R', 'R0 D Z mu alpha')
    generate_country_model('seiir_int_reparam', 'S E Ir Iu R', 'R0 D Z mu alpha tact dtact kbeta')
    generate_country_model('seiir_intexp_reparam', 'S E Ir Iu R', 'R0 D Z mu alpha tact k')
    generate_country_model('seiir_intrem_reparam', 'S E Ir Iu R', 'R0 D Z mu alpha tact l')

    # add model here as argument
    generate_country('sir', 'sir_int', 'sir_int_r0', 'seir', 'seir_int', 'seiir', 'seiir_int',
                     'sir_reparam', 'sir_int_reparam', 'sir_intexp_reparam', 'sir_intrem_reparam', 'seir_reparam', 'seir_int_reparam', 'seir_intexp_reparam', 'seir_intrem_reparam', 'seiir_reparam', 'seiir_int_reparam','seiir_intexp_reparam', 'seiir_intrem_reparam')
    generate_canton('sei_c', 'seii_c', 'seiin', 'seiin_interventions')


if __name__ == '__main__':
    main()
