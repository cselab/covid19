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
f_lambda = dict()
f_lambda_int = dict()
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
    f_lambda[folder] = beta - gamma
    f_lambda_int[folder] = beta_int - gamma
    f_infected[folder] = js['y-data'][-1]


def Col(dictx, dicty):
    x = []
    y = []
    for k in dictx:
        if k in dicty:
            x.append(dictx[k])
            y.append(dicty[k])
    return x, y


'''
fig, [ax, ax2] = plt.subplots(2)
ax.set_title(r"$\lambda*=\beta - \gamma$ before intervention")
ax.scatter(*Col(f_lambda, f_infected))
ax.axvline(x=0, color='black', linestyle=':')
ax.set_ylabel("infected")

ax2.set_title(r"$\lambda*=\beta - \gamma$ after intervention")
ax2.scatter(*Col(f_lambda_int, f_infected))
ax2.axvline(x=0, color='black', linestyle=':')
ax2.set_ylabel("infected")

fig.tight_layout()
fig.savefig("scatter.pdf")
'''

displayname = {
    "RussianFederation": "Russia",
    "UnitedKingdom": "UK",
    "BosniaandHerzegovina": "BiH",
    "NorthMacedonia": "North Macedonia",
}

fig, axes = plt.subplots(1, figsize=(9, 6))
axes = [axes, axes]

color = dict()

for f, ax in zip([f_R0, f_R0_int], axes):
    before = (f == f_R0)
    ax.axvline(x=1, color='black', linestyle=':')
    i = 0
    cc = np.array(list(f_R0.keys()))
    vv = np.array(list(f_R0.values()))
    arg = np.argsort(vv)
    ax.get_yaxis().set_visible(False)
    for i, a in enumerate(arg):
        c = cc[a]
        xy = f[c], i
        p = ax.scatter(*xy, s=16, c=color.get(c, None))
        color[c] = p.get_facecolor()
        ax.annotate(displayname.get(c, c),
                    xy=xy,
                    fontsize=5,
                    xytext=(4, 0),
                    textcoords='offset points',
                    va='center')
    ax.set_xlim(0, 3)
    ax.text(0.03 if before else 0.55,
            1.01,
            r"$R_0$ {:} intervention".format("before" if before else "after"),
            transform=ax.transAxes, fontsize=15)
fig.tight_layout()
fig.savefig("scatter_R0.pdf")
