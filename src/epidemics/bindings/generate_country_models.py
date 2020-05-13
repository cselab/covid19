#!/usr/bin/env python3

import jinja2
import os
PATH = os.path.normpath(os.path.dirname(__file__))

NO_EDIT_NOTE = f"""
// THIS FILE WAS CREATED FROM {os.path.basename(__file__)}
// DO NOT EDIT!

""".lstrip()

with open(os.path.join(PATH, 'country_model.template.cpp')) as f:
    COUNTRY_TEMPLATE = jinja2.Template(f.read())


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


def generate_country_model(name, state, params):
    """Generate bindings for one country model."""
    state = state.split()
    params = params.split()
    code = NO_EDIT_NOTE + COUNTRY_TEMPLATE.render(
            NAME=name, STATE=state, PARAMS=params)
    save_if_modified(os.path.join(PATH, f'_country.{name}.generated.cpp'), code)


def generate_country(*models):
    """Generate _country.generated.cpp, responsible for libepidemics.country submodule."""
    with open(os.path.join(PATH, 'country.template.cpp')) as f:
        template = jinja2.Template(f.read())
    code = NO_EDIT_NOTE + template.render(MODELS=models)
    save_if_modified(os.path.join(PATH, '_country.generated.cpp'), code)


def main():
    generate_country_model('sir', 'S I R', 'beta gamma')
    generate_country_model('seiir', 'S E Ir Iu R', 'beta mu alpha Z D')
    generate_country('sir', 'seiir')


if __name__ == '__main__':
    main()
