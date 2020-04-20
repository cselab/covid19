#!/usr/bin/env python3
# Author: Pantelis R. Vlachas
# Date:   20/4/2020
# Email:  pvlachas@ethz.ch
from pathlib import Path
import numpy as np

import json
import os
import pickle
import sys


CURRENT_DIR = Path(__file__)
# print(CURRENT_DIR)
# print(os.path.dirname(__file__))
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# print(CURRENT_DIR)
PROJECT_DIR = Path(CURRENT_DIR).parent
PROJECT_DIR = PROJECT_DIR.parent
print(PROJECT_DIR)
sys.path.append(PROJECT_DIR) 
# DATA_DIR = Path(CURRENT_DIR).parent / 'data'
# print(DATA_DIR)
# TOOLS_DIR = Path(CURRENT_DIR).parent / 'tools'
# print(TOOLS_DIR)
# sys.path.insert(0, TOOLS_DIR) 


# CURRENT_DIR = Path(__file__)
# print(CURRENT_DIR)
# print(os.path.dirname(__file__))
# CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# print(CURRENT_DIR)
# PROJECT_DIR = Path(CURRENT_DIR).parent
# print(PROJECT_DIR)
# DATA_DIR = Path(CURRENT_DIR).parent / 'data'
# print(DATA_DIR)
# sys.path.insert(0, PROJECT_DIR) 
# TOOLS_DIR = Path(CURRENT_DIR).parent / 'tools'
# print(TOOLS_DIR)
# sys.path.insert(0, TOOLS_DIR) 

# from tools.io import *
# from io import *
from epidemics.tools.io import download_and_save
# import epidemics.epidemics

# DATA_DOWNLOADS_DIR = DATA_DIR / 'downloads'
# print(DATA_DOWNLOADS_DIR)


# def fetch_openzh_covid_data(*, cache_duration=3600):
#     """
#     Returns a dictionary of lists {canton abbreviation: number of cases per day}.
#     """
#     url = 'https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv'
#     path = DATA_DOWNLOADS_DIR / 'covid19_cases_switzerland_openzh.csv'

#     raw = download_and_save(url, path, cache_duration=cache_duration)
#     rows = raw.decode('utf8').split()
#     cantons = rows[0].split(',')[1:-1]  # Skip the "Date" and "CH" cell.

#     data = {canton: [] for canton in cantons}
#     for day in rows[1:]:  # Skip the header.
#         cells = day.split(',')[1:-1]  # Skip "Date" and "CH".
#         assert len(cells) == len(cantons), (len(cells), len(cantons))

#         for canton, cell in zip(cantons, cells):
#             data[canton].append(float(cell or 'nan'))
#     return data


# cases_per_country = fetch_openzh_covid_data()

# print(cases_per_country)

# # print(ark)






