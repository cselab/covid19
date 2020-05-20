#!/usr/bin/env python3

import json
import numpy as np
import argparse
import os
from datetime import datetime
from glob import glob
import re
import pandas as pd

# XXX path to folder with output from `request_country.py`
datafolder = "."

LAST_DAY = "2020-05-18"

def get_folder(country):
    return os.path.join(datafolder, country.replace(' ', ''))

def get_foldername(country):
    return country.replace(' ', '')


def folder_to_country(folder):
    displayname = {
        "RussianFederation": "Russia",
        "UnitedKingdom": "United Kingdom",
        "BosniaandHerzegovina": "Bosnia and Herzegovina",
        "NorthMacedonia": "North Macedonia",
        "CzechRepublic": "Czechia",
    }
    return displayname.get(folder, folder)

# folder name to parameter
f_tact = dict()
f_tact_day = dict()
f_day0 = dict()   # day of the first row of `x-data`

def days_to_delta(t):
    return np.timedelta64(int(t + 0.5), 'D')

folders = []
for path in glob(os.path.join(datafolder, "*", "intervals.json")):
    folder = re.findall(os.path.join(datafolder, "(.*)", "intervals.json"),
                        path)[0]
    with open(path) as f:
        js = json.loads(f.read())
    p = js['mean_params']
    gamma = p['gamma']
    kbeta = p['kbeta']
    beta = p['R0'] * gamma
    beta_int = beta * kbeta
    folders.append(folder)
    day0 = np.datetime64(LAST_DAY) - days_to_delta(max(js['x-data']))
    f_tact[folder] = p['tact']
    f_day0[folder] = day0
    f_tact_day[folder] = day0 + days_to_delta(p['tact'])


for c in sorted(folders):
    ts = pd.to_datetime(str(f_tact_day[c]))
    #day = ts.strftime('%Y.%m.%d')
    day = ts.strftime('%B %d')
    print("{:},{:}".format(folder_to_country(c), day))
