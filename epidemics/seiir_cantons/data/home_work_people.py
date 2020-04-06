#!/usr/bin/env python3

import numpy as np
import itertools
import json
import collections

fname = 'commute.csv'
d = np.genfromtxt(fname, delimiter=',', names=True, dtype=None, encoding=None)

home = d['home']
work = d['work']
people = d['people']

codes = np.unique(home)

matrix = collections.defaultdict(dict)

for c0,c1 in itertools.product(codes, repeat=2):
    matrix[c0][c1] = int(people[np.where((home == c0) & (work == c1))].sum())

oname = 'home_work_people.json'
with open(oname, 'w') as o:
    json.dump(matrix, o)

