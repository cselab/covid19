#!/usr/bin/env python

import matplotlib.pyplot as plt
import json
import os
from glob import glob
import re
import numpy as np

# XXX path to folder with output from `request_country.py`
datafolder = "."


def get_folder(country):
    return os.path.join(datafolder, country.replace(' ', ''))


def get_foldername(country):
    return country.replace(' ', '')


# folder name to parameter
f_R0 = dict()
f_R0_int = dict()
f_tact = dict()
f_infected = dict()

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
    f_R0[folder] = p['R0']
    f_R0_int[folder] = p['R0'] * kbeta
    f_tact[folder] = p['tact']
    f_infected[folder] = js['y-data'][-1]


def Col(dictx, dicty):
    x = []
    y = []
    for k in dictx:
        if k in dicty:
            x.append(dictx[k])
            y.append(dicty[k])
    return x, y


displayname = {
    "RussianFederation": "Russia",
    "UnitedKingdom": "UK",
    "BosniaandHerzegovina": "BiH",
    "NorthMacedonia": "North Macedonia",
    "CzechRepublic": "Czechia",
}

fig, axes = plt.subplots(2, figsize=(9, 6))

color = dict()

for f, ax in zip([f_R0, f_R0_int], axes):
    before = (f == f_R0)
    ax.axvline(x=1, color='black', linestyle=':', zorder=-10)
    i = 0
    for c in f:
        xy = f[c], f_tact[c]
        p = ax.scatter(*xy, s=16, c=color.get(c, None))
        color[c] = p.get_facecolor()
        ax.annotate(displayname.get(c, c),
                    xy=xy,
                    fontsize=7,
                    xytext=(4, 0),
                    textcoords='offset points',
                    va='center')
    ax.set_xlim(-0.05, 3.05)
    ax.text(0.03 if not before else 0.55,
            1.01,
            r"$R_0$ {:} intervention".format("before" if before else "after"),
            transform=ax.transAxes, fontsize=15)
    ax.set_ylabel('intervention time')
fig.tight_layout()
fig.savefig("scatter_tact.pdf")
