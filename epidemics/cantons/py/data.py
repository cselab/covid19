#!/usr/bin/env python3

from datetime import datetime
import json
import math
import numpy as np
import os
import urllib.request

DATA_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data'))
MIJ_MATRIX_JSON = os.path.join(DATA_DIR, 'home_work_people.json')

CANTON_POPULATION = dict(zip(
    'ZH BE LU UR SZ OW NW GL ZG FR SO BS BL SH AR AI SG GR AG TG TI VD VS NE GE JU'.split(),
    [
        1520968,
        1034977,
        409557,
        36433,
        159165,
        37841,
        43223,
        40403,
        126837,
        318714,
        273194,
        200298,
        289527,
        81991,
        55234,
        16145,
        507697,
        198379 ,
        678207,
        276472,
        353343,
        799145,
        343955,
        176850,
        499480,
        73419,
    ]))


def fetch(cache=True):
    fname = os.path.join(DATA_DIR, 'covid19_cases_switzerland_openzh.csv')
    if not cache or not os.path.isfile(fname):
        url = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv"
        req = urllib.request.urlopen(url)
        data = req.read()
        with open(fname, 'wb') as f:
            f.write(data)
        return data.decode('utf8')
    with open(fname) as f:
        return f.read()


def fetch_canton_data(cache=True):
    """Fetch per-day per-canton number of cases.

    Returns a tuple of two values:
        - canton map {abbreviation string: index}
        - matrix m[day][canton index], number of cases (list of lists of floats)
    """
    rows = fetch(cache=cache).split()
    cantons = rows[0].split(',')[1:-1]  # Skip the "Date" cell and "CH",


    A = set(cantons)
    B = set(CANTON_POPULATION.keys())
    assert A == B, (A - B, B - A)

    data = []
    for day in rows[1:]:
        day = day.split(',')[1:-1]  # Skip "Date" and "CH".
        assert len(day) == len(cantons), (len(day), len(cantons))

        data.append([float(cell or 'nan') for cell in day])

    # data = np.array(data)

    return {key: k for k, key in enumerate(cantons)}, data


def get_raw_Mij(cantons):
    """Returns a matrix M[i][j] which represents daily commute from i to j.

    Arguments:
        cantons: List of canton abbreviations.
    """
    with open(MIJ_MATRIX_JSON) as f:
        Mjson = json.load(f)
    for canton in cantons:
        print(canton, Mjson[canton][canton], CANTON_POPULATION[canton])
        Mjson[canton][canton] = 0
    Mij = [[Mjson[a][b] for a in cantons] for b in cantons]

    return Mij


def get_symmetric_Mij(cantons):
    N = len(cantons)
    Mij = get_raw_Mij(cantons)
    Mij = [[Mij[i][j] + Mij[j][i] for j in range(N)] for i in range(N)]
    return Mij

def prepare_data_for_cpp():
    """Generate ../data/cantons_data.dat, used by CantonsData class in the C++ code.

    File format:
        <number of cantons N>
        abbreviation1 ... abbreviationN
        population1 ... populationN

        M_11 ... M_1N
        ...
        M_N1 ... M_NN

        <number of known data points M>
        day1 canton_index1 number_of_infected1
        ...
        dayM canton_indexM number_of_infectedM

    The known data points are all known values of numbers of infected.
    Note that covid19_cases_switzerland_openzh.csv (see `fetch`) has many missing values.
    """

    canton_to_index, infected = fetch_canton_data()
    Mij = get_raw_Mij(canton_to_index)

    with open(os.path.join(DATA_DIR, 'cantons_data.dat'), 'w') as f:
        f.write(str(len(canton_to_index)) + '\n')
        f.write(" ".join(canton_to_index) + '\n')
        f.write(" ".join(str(CANTON_POPULATION[c]) for c in canton_to_index) + '\n\n')

        for row in Mij:
            f.write(" ".join(str(x) for x in row) + '\n')
        f.write('\n')

        known_data_points = [
            (d, c, country_day_value)
            for d, day_values in enumerate(infected)
            for c, country_day_value in enumerate(day_values)
            if not math.isnan(country_day_value)
        ]
        f.write(str(len(known_data_points)) + '\n')
        for data_point in known_data_points:
            f.write("{} {} {}\n".format(*data_point))


if __name__ == '__main__':
    prepare_data_for_cpp()
