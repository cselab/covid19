#!/usr/bin/env python3

import numpy as np

from datetime import datetime
import json
import math
import os
import sys
import urllib.request

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'build'))

try:
    import libsolver
except ModuleNotFoundError:
    sys.exit("libsolver not found. Did you forget to compile the C++ code?")

from epidemics.data import DATA_CACHE_DIR
from epidemics.data.cases import get_region_cases
from epidemics.tools.tools import flatten
import epidemics.data.swiss_cantons as swiss_cantons


class ModelData:
    """Model data such as canton population and Mij matrix.

    For conveniece, we store parameters beta, mu, alpha and other scalars
    separately (as libsolver.Parameters).

    Arguments:
        canton_keys: List of canton abbreviated names
        canton_population: List of population size of corresponding cantons.
        Mij: A numpy matrix of canton-canton number of commuters.
        external_cases: A matrix [day][canton] of estimated number of foreign
                        infected people visiting given canton at given day.
    """
    def __init__(self, canton_keys, canton_population, Mij, external_cases):
        self.num_cantons = len(canton_keys)
        self.canton_keys = canton_keys
        self.canton_population = canton_population
        self.Mij = Mij
        self.external_cases = external_cases

    def to_cpp(self):
        """Return the libsolver.ModelData instance.

        Needed when running the model from Python using the C++ implementation."""
        keys = {key: k for k, key in enumerate(self.canton_keys)}
        return libsolver.ModelData(self.num_cantons, keys, self.canton_population,
                                   flatten(self.Mij), flatten(self.external_cases))

    def save_cpp_dat(self, path=DATA_CACHE_DIR / 'cpp_model_data.dat'):
        """Generate cpp_model_data.dat, the data for the C++ ModelData class.

        Needed when running Korali from C++, when `to_cpp` is not available.

        File format:
            <number of cantons N>
            abbreviation1 ... abbreviationN
            population1 ... populationN

            M_11 ... M_1N
            ...
            M_N1 ... M_NN

            <number of days D for external cases>
            <external cases for day 1 for canton 1> ... <external cases for day 1 for canton N>
            ...
            <external cases for day D for canton 1> ... <external cases for day D for canton N>
        """
        with open(path, 'w') as f:
            f.write(str(self.num_cantons) + '\n')
            f.write(' '.join(self.canton_keys) + '\n')
            f.write(' '.join(str(p) for p in self.canton_population) + '\n\n')

            for row in self.Mij:
                f.write(' '.join(str(x) for x in row) + '\n')
            f.write('\n')

            f.write(str(len(self.external_cases)) + '\n')
            for day in self.external_cases:
                f.write(' '.join(str(x) for x in day) + '\n')
        print(f"Stored model data to {path}.")


class ValidationData:
    """Measured data that the model is predicting."""
    def __init__(self, canton_keys, cases_per_country):
        self.canton_keys = canton_keys
        self.cases_per_country = cases_per_country
        self.canton_to_index = {key: k for k, key in enumerate(canton_keys)}

        # Not all elements of the `cases_per_country` matrix are known, so we
        # create a list of tuples (day, canton index, number of cases).
        self.cases_data_points = [
            (d, self.canton_to_index[c], day_value)
            for c, canton_values in cases_per_country.items()
            for d, day_value in enumerate(canton_values)
            if not math.isnan(day_value)
        ]

    def save_cpp_dat(self, path=DATA_CACHE_DIR / 'cpp_validation_data.dat'):
        """Generate cpp_validation_data.dat, the data for the C++ ValidationData class.

        File format:
            <number of known data points M>
            day1 canton_index1 number_of_cases1
            ...
            dayM canton_indexM number_of_casesM

        The known data points are all known values of numbers of cases. The
        canton index refers to the order in cpp_model_data.dat
        Note that covid19_cases_switzerland_openzh.csv (see `fetch`) has many missing values.
        """
        with open(path, 'w') as f:
            f.write(str(len(self.cases_data_points)) + '\n')
            for data_point in self.cases_data_points:
                f.write('{} {} {}\n'.format(*data_point))
        print(f"Stored validation data to {path}.")


def get_model_data(include_foreign=True):
    """Creates the ModelData instance with default data."""
    keys = swiss_cantons.CANTON_KEYS_ALPHABETICAL
    population = [swiss_cantons.CANTON_POPULATION[c] for c in keys]

    Mij = swiss_cantons.get_default_Mij_numpy(keys)

    if include_foreign:
        swiss_cases = get_region_cases('switzerland')
        external_cases = swiss_cantons.get_external_cases(
                swiss_cases.get_date_of_first_confirmed(),
                num_days=len(swiss_cases.confirmed) + 10)
        # A matrix [d][c] of foreign infected people visiting canton c at day d.
        external_cases = [external_cases[c] for c in keys]
    else:
        external_cases = [[] for c in keys]

    return ModelData(keys, population, Mij, external_cases)


def get_validation_data():
    keys = swiss_cantons.CANTON_KEYS_ALPHABETICAL
    cases_per_country = swiss_cantons.fetch_openzh_covid_data()
    return ValidationData(keys, cases_per_country)


if __name__ == '__main__':
    make_model_data().save_cpp_dat()
    make_validation_data().save_cpp_dat()
