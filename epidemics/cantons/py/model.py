#!/usr/bin/env python3

import numpy as np

from datetime import datetime
import json
import math
import os
import sys
import urllib.request

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'build'))

try:
    import libepidemics
except ModuleNotFoundError:
    sys.exit("libepidemics not found. Did you forget to compile the C++ code?")

from epidemics.data import DATA_CACHE_DIR
from epidemics.data.cases import get_region_cases
from epidemics.tools.tools import flatten
import epidemics.data.swiss_cantons as swiss_cantons
import epidemics.data.swiss_municipalities as swiss_municipalities


class PyDesignParameters:
    """Design parameters such as region population and Mij matrix.

    To be used from Python and finally converted to the C++ DesignParameters
    when creating a C++ Solver.

    Arguments:
        region_keys: List of region names.
        region_population: List of population size of corresponding regions.
        Mij: A numpy matrix of region-region number of commuters.
        ext_com_Iu: A matrix [day][region] of estimated number of foreign
                    infected people visiting given region at given day.
        Ui: User-defined, shape (K)
    """
    def __init__(self, region_keys, region_population, Mij, Cij, *, ext_com_Iu=[], Ui=[]):
        K = len(region_keys)
        assert len(region_population) == K
        assert Mij.shape == (K, K)
        assert Cij.shape == (K, K)
        assert all(len(row) == K for row in ext_com_Iu), \
                (K, [len(row) for row in ext_com_Iu])
        if not len(Ui):
            Ui = [0] * K
        assert len(Ui) == K

        self.num_regions = K
        self.region_keys = region_keys
        self.region_population = region_population
        self.Mij = Mij
        self.Cij = Cij
        self.ext_com_Iu = ext_com_Iu
        self.Ui = Ui

        self.key_to_index = {key: k for k, key in enumerate(region_keys)}

    def to_cpp(self):
        """Return the libepidemics.DesignParameters instance.

        Needed when running the model from Python using the C++ implementation."""
        return libepidemics.cantons.DesignParameters(
                self.region_keys, self.region_population,
                flatten(self.Mij), flatten(self.Cij),
                flatten(self.ext_com_Iu), self.Ui)

    def save_cpp_dat(self, path=DATA_CACHE_DIR / 'cpp_design_parameters.dat'):
        """Generate cpp_design_parameters.dat, the data for the C++ DesignParameters class.

        Needed when running Korali from C++, when `to_cpp` is not available.

        File format:
            <number of regions N>
            abbreviation1 ... abbreviationN
            population1 ... populationN

            M_11 ... M_1N
            ...
            M_N1 ... M_NN

            <number of days D for external cases>
            <external cases for day 1 for region 1> ... <external cases for day 1 for region N>
            ...
            <external cases for day D for region 1> ... <external cases for day D for region N>
        """
        with open(path, 'w') as f:
            f.write(str(self.num_regions) + '\n')
            f.write(' '.join(self.region_keys) + '\n')
            f.write(' '.join(str(p) for p in self.region_population) + '\n\n')

            for row in self.Mij:
                f.write(' '.join(str(x) for x in row) + '\n')
            f.write('\n')

            f.write(str(len(self.ext_com_Iu)) + '\n')
            for day in self.ext_com_Iu:
                f.write(' '.join(str(x) for x in day) + '\n')
            f.write(' '.join(str(u) for u in self.Ui) + '\n')
        print(f"Stored design parameters to {path}.")


class ReferenceData:
    """Measured data that the model is predicting."""
    def __init__(self, region_keys, cases_per_country):
        self.region_keys = region_keys
        self.cases_per_country = cases_per_country
        self.region_to_index = {key: k for k, key in enumerate(region_keys)}

        # Not all elements of the `cases_per_country` matrix are known, so we
        # create a list of tuples (day, region index, number of cases).
        self.cases_data_points = [
            (d, self.region_to_index[c], day_value)
            for c, region_values in cases_per_country.items()
            for d, day_value in enumerate(region_values)
            if not math.isnan(day_value)
        ]

    def save_cpp_dat(self, path=DATA_CACHE_DIR / 'cpp_reference_data.dat'):
        """Generate cpp_reference_data.dat, the data for the C++ ReferenceData class.

        File format:
            <number of known data points M>
            day1 region_index1 number_of_cases1
            ...
            dayM region_indexM number_of_casesM

        The known data points are all known values of numbers of cases. The
        region index refers to the order in cpp_design_parameters.dat
        Note that covid19_cases_switzerland_openzh.csv (see `fetch`) has many missing values.
        """
        with open(path, 'w') as f:
            f.write(str(len(self.cases_data_points)) + '\n')
            for data_point in self.cases_data_points:
                f.write('{} {} {}\n'.format(*data_point))
        print(f"Stored reference data to {path}.")


def get_cantom_design_parameters(include_foreign=True):
    """Creates the PyDesignParameters instance with default data."""
    keys = swiss_cantons.CANTON_KEYS_ALPHABETICAL
    population = [swiss_cantons.CANTON_POPULATION[c] for c in keys]

    Mij = swiss_cantons.get_Mij_numpy(keys)
    Cij = swiss_cantons.get_Cij_numpy(keys)

    if include_foreign:
        swiss_cases = get_region_cases('switzerland')
        num_days = len(swiss_cases.confirmed) + 10
        ext_com_Iu = swiss_cantons.get_external_Iu(
                swiss_cases.get_date_of_first_confirmed(), num_days=num_days)
        # Transpose from {canton: [day1, ...]} to [day][canton].
        ext_com_Iu = [[ext_com_Iu[c][d] for c in keys] for d in range(num_days)]
    else:
        ext_com_Iu = []  # Data for 0 days.

    return PyDesignParameters(keys, population, Mij, Cij, ext_com_Iu=ext_com_Iu)


def get_canton_model_data(*args, **kwargs):
    print("WARNING: get_canton_model_data was renamed to get_canton_design_parameters!")
    return get_canton_design_parameters(*args, **kwargs)


def get_canton_reference_data():
    keys = swiss_cantons.CANTON_KEYS_ALPHABETICAL
    cases_per_country = swiss_cantons.fetch_openzh_covid_data()
    return ReferenceData(keys, cases_per_country)


def get_municipality_design_parameters():
    namepop = swiss_municipalities.get_name_and_population()
    commute = swiss_municipalities.get_commute()

    # TODO: Comparison with reference data requires aggregation wrt cantons,
    # since that's the only reference data we have.
    # cantons = swiss_municipalities.get_municipality_commute()

    key_to_index = {key: k for k, key in enumerate(namepop['key'])}
    N = len(key_to_index)
    Cij = np.zeros((N, N))
    for key_home, key_work, num_people in zip(
            commute['key_home'],
            commute['key_work'],
            commute['num_people']):
        home = key_to_index.get(key_home)
        work = key_to_index.get(key_work)
        if home is None or work is None:
            continue
        Cij[work, home] += num_people

    # NOTE: This Mij is wrong.
    Mij = Cij + Cij.transpose()

    return PyDesignParameters(namepop['key'], namepop['population'], Mij, Cij)


def get_municipality_model_data():
    print("WARNING: get_municipality_model_data was renamed to get_municipality_design_parameters!")
    return get_municipality_design_parameters()


if __name__ == '__main__':
    get_canton_design_parameters().save_cpp_dat()
    get_canton_reference_data().save_cpp_dat()
    # get_municipality_design_parameters().save_cpp_dat()
