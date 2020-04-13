#!/usr/bin/env python3

import numpy as np

import collections
import itertools
import json
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data import DATA_CACHE_DIR

fname = DATA_CACHE_DIR / 'switzerland_commute_admin_ch.csv'
d = np.genfromtxt(fname, delimiter=',', names=True, dtype=None, encoding=None)

home = d['home']
work = d['work']
people = d['people']

codes = np.unique(home)

matrix = collections.defaultdict(dict)

for c0,c1 in itertools.product(codes, repeat=2):
    matrix[c0][c1] = int(people[np.where((home == c0) & (work == c1))].sum())

# TODO: Move to epidemics/data/swiss_cantons.py!
oname = 'home_work_people.json'
with open(oname, 'w') as o:
    json.dump(matrix, o)

