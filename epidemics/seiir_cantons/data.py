import numpy as np
import os
import urllib.request
from datetime import datetime


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
    fname = "covid19_cases_switzerland_openzh.csv"
    if not cache or not os.path.isfile(fname):
        url = "https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv"
        req = urllib.request.urlopen(url)
        data = req.read()
        with open(fname, 'wb') as f:
            f.write(data)
        return data
    with open(fname) as f:
        return f.read()


def fetch_canton_data(cache=True):
    """Fetch per-day per-canton number of cases.

    Returns a tuple of two values:
        - canton abbreviations (list of strings)
        - matrix m[day][canton index], number of cases (list of lists of floats)
    """
    rows = fetch(cache=cache).split()
    cantons = rows[0].split(',')[1:]  # Skip the "Date" cell.

    assert set(cantons) == set(CANTON_POPULATION.keys()), \
            set(cantons) ^ set(CANTON_POPULATION.keys())

    data = []
    for day in rows[1:]:
        day = day.split(',')[1:]  # Skip the date.
        assert len(day) == len(cantons), (len(day), len(cantons))

        data.append([float(cell or 'nan') for cell in day])
        print(data[-1])

    # data = np.array(data)

    return cantons, data


if __name__ == '__main__':
    fetch_canton_data()
