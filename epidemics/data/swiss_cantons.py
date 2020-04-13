"""
This file provides access to all relevant data aobut Swiss cantons and COVID-19:
    - list of cantons and their population count
    - connectivity matrix between cantons (number of daily commuters)
    - number of cases per day, for each canton


"""

from epidemics.data import DATA_CACHE_DIR, DATA_DOWNLOADS_DIR, DATA_FILES_DIR
from epidemics.tools.cache import cache_to_file
from epidemics.tools.io import download_and_save
import numpy as np

# https://en.wikipedia.org/wiki/Cantons_of_Switzerland
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

CANTON_KEYS_ALPHABETICAL = sorted(CANTON_POPULATION.keys())

def fetch_openzh_covid_data(*, cache_duration=3600):
    """
    Returns a dictionary of lists {canton abbreviation: number of cases per day}.
    """
    url = 'https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv'
    path = DATA_DOWNLOADS_DIR / 'covid19_cases_switzerland_openzh.csv'

    raw = download_and_save(url, path, cache_duration=cache_duration)
    rows = raw.decode('utf8').split()
    cantons = rows[0].split(',')[1:-1]  # Skip the "Date" and "CH" cell.

    data = {canton: [] for canton in cantons}
    for day in rows[1:]:  # Skip the header.
        cells = day.split(',')[1:-1]  # Skip "Date" and "CH".
        assert len(cells) == len(cantons), (len(cells), len(cantons))

        for canton, cell in zip(cantons, cells):
            data[canton].append(float(cell or 'nan'))
    return data



COMMUTE_ADMIN_CH_CSV = DATA_FILES_DIR / 'switzerland_commute_admin_ch.csv'

@cache_to_file(DATA_CACHE_DIR / 'home_work_people.json',
               dependencies=[COMMUTE_ADMIN_CH_CSV])
def get_Mij_home_work_admin_ch_json():
    """
    Returns a dictionary
    {canton1: {canton2: number of commuters between canton1 and canton2, ...}, ...}.
    """
    cantons = set()
    entries = []
    with open(COMMUTE_ADMIN_CH_CSV) as f:
        header = f.readline()
        for line in f:
            home, work, people = line.split(',')
            if work == 'ZZ':
                continue
            cantons.add(home)
            cantons.add(work)
            entries.append((home, work, int(people)))

    Mij = {c1: {c2: 0 for c2 in cantons} for c1 in cantons}
    for home, work, people in entries:
        Mij[home][work] += people
        Mij[work][home] += people

    return Mij


def Mij_json_to_numpy(json, canton_order):
    """Returns the number of commuters data as a matrix.

    Arguments:
        json: A commute matrix in a JSON dictionary format (see get_Mij_home_work_admin_ch_json).
        canton_order: The desired row and column order in the output matrix.
    """
    assert len(canton_order) == len(json)
    out = np.zeros((len(canton_order), len(canton_order)))
    for index1, c1 in enumerate(canton_order):
        for index2, c2 in enumerate(canton_order):
            out[index1][index2] = json[c1][c2]
    return out


get_default_Mij_json = get_Mij_home_work_admin_ch_json

def get_default_Mij_numpy(canton_order):
    """Return the Mij numpy matrix using the data from the default source (admin.ch)."""
    return Mij_json_to_numpy(get_default_Mij_json(), canton_order)
