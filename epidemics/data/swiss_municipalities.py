import numpy as np
import pandas as pd

import re
import os

from epidemics.data import DATA_CACHE_DIR, DATA_DOWNLOADS_DIR
from epidemics.tools.cache import cache, cache_to_file
from epidemics.tools.io import download_and_save

# Note of code design:
# 1) "Get" functions are split into smaller ones, each on of takes cares of
#    certain "column" of data. That way we abstract away the data source,
#    which simplifies how we replace sources with another ones.
#
# 2) For example, the function `get_municipality_cantons` may use a source
#    that is different from the source for `Mij`.
#
# 3) Cached functions are those that take a bunch of time
#    (those that need access to Excel files).
#
# Terminology:
#     region: A general term for countries / cantons / municipalities.
#     municipality number: Official 4-digit municipality number, sometimes 0-padded, sometimes not.
#     key: A unique ID of a region. For municipalities we use 'MUN-0123', where 0123 is the 0-padded municipality number.

@cache
def bfs_residence_work_xls(header=(3, 4), usecols=None):
    """Return the residence-workplace commute Excel file as a pandas.DataFrame."""
    url = 'https://www.bfs.admin.ch/bfsstatic/dam/assets/8507281/master'
    path = DATA_DOWNLOADS_DIR / 'bfs_residence_work.xlsx'
    download_and_save(url, path, load=False)
    print(f"Loading {path}...", flush=True)
    sheet = pd.read_excel(path, sheet_name='Commune of residence perspect.',
                          header=header, skipfooter=4, usecols=usecols)
    return sheet


@cache
def bfs_population_xls(usecols=None):
    """Return the municipality Excel file as a pandas.DataFrame."""
    url = 'https://www.bfs.admin.ch/bfsstatic/dam/assets/9635941/master'
    path = DATA_DOWNLOADS_DIR / 'bfs_municipality_population.xlsx'
    download_and_save(url, path, load=False)
    print(f"Loading {path}...", flush=True)
    sheet = pd.read_excel(path, sheet_name='2018',
                          header=1, skipfooter=5, usecols=usecols)
    return sheet


@cache
@cache_to_file(DATA_CACHE_DIR / 'bfs_residence_work_cols12568.df.csv')
def get_residence_work_cols12568():
    # (residence canton initial,
    #  residence commune number,
    #  work canton initial,
    #  work commune number,
    #  number of employed people)
    df = bfs_residence_work_xls(header=4, usecols=(1, 2, 5, 6, 8))
    df.columns = ('canton_home', 'number_home', 'canton_work', 'number_work', 'num_people')
    return df


def _number_to_key(number):
    """
    Add a prefix 'MUN-' to avoid confusion between the
    municipality numbers and indices. The numbers go
    from 1 to ~6800, while indices from 0 to ~2200.
    """
    return 'MUN-{:04}'.format(int(number))


@cache_to_file(DATA_CACHE_DIR / 'bfs_municipality_namepop.df.csv')
def get_name_and_population():
    """Returns a pandas DataFrame with municipality names and population.

    >>> get_municipality_names_and_population()
               key                name  population
    0     MUN-0001     Aeugst am Albis        1982
    1     MUN-0002  Affoltern am Albis       12229
    ...        ...                 ...         ...
    """
    sheet = bfs_population_xls(usecols=('Region', 'Total'))

    rows = []

    # We are matching strings of the norm '......1234 Municipality Name'.
    NUMBER_NAME_RE = re.compile(r'\.+(\d+) (.+)')
    for region, total in zip(sheet['Region'], sheet['Total']):
        match = NUMBER_NAME_RE.match(region)
        if not match:
            continue  # Not a municipality, but an aggregate.

        key = _number_to_key(match.group(1))
        name = match.group(2)
        rows.append((key, name, total))

    return pd.DataFrame(rows, columns=('key', 'name', 'population'))


def get_cantons():
    """Returns a DataFrame with municipality cantons.

    >>> get_cantons()
               key canton
    0     MUN-0001     ZH
    1     MUN-0002     ZH
    ...        ...    ...
    """
    commute = get_residence_work_cols12568()
    print(commute)
    return pd.DataFrame({
        'key': list(map(_number_to_key, commute['number_home'])),
        'canton': commute['canton_home'],
    })


def get_commute():
    """Returns a DataFrame with data on commute between municipality.

    >>> get_commute()
           key_home  key_work  num_people
    0      MUN-0001  MUN-0001         147
    1      MUN-0001  MUN-0002         106
    ...         ...       ...         ...
    """
    commute = get_residence_work_cols12568()
    return pd.DataFrame({
        'key_home': list(map(_number_to_key, commute['number_home'])),
        'key_work': list(map(_number_to_key, commute['number_work'])),
        'num_people': commute['num_people']
    })

def get_shape_file():
    """Downloads and returns path to shape file with multicipalities.

    >>> get_shape_file()
    DIR/swissBOUNDARIES3D_1_3_TLM_BEZIRKSGEBIET
    """
    from zipfile import ZipFile
    zippath = DATA_DOWNLOADS_DIR / "swissBOUNDARIES3D.zip"
    download_and_save("https://shop.swisstopo.admin.ch/shop-server/resources/products/swissBOUNDARIES3D/download", zippath)

    DATA_MAP_DIR = DATA_DOWNLOADS_DIR / "map"
    os.makedirs(DATA_MAP_DIR, exist_ok=True)

    shapefile = "BOUNDARIES_2020/DATEN/swissBOUNDARIES3D/SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_3_TLM_BEZIRKSGEBIET"
    with ZipFile(zippath, 'r') as zipobj:
        for member in zipobj.namelist():
            if shapefile in member:
                basename = member
                path = DATA_MAP_DIR / os.path.basename(member)
                if not os.path.isfile(path):
                    print("extracting '{:}'".format(basename))
                    with open(path, 'wb') as f:
                        f.write(zipobj.read(member))
    return os.path.splitext(path)[0]

