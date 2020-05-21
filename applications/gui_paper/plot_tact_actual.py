#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import json
import os
from glob import glob
import re
import numpy as np
import pandas as pd

# XXX path to folder with output from `request_country.py`
datafolder = "."

LAST_DAY = "2020-05-18"

csv = "data/delay.csv"
csv = pd.read_csv(csv)

print(len(csv))
"""
for i in range(len(out)):
  date = merged['date'][i]
  exact = merged['exactDate'][i]
  if exact and exact != 'nan':
      date = pd.to_datetime(date)
      exact = pd.to_datetime(exact)
      delay = date - exact
      if not pd.isna(delay):
          out['delay'][i] = "{:+d}".format(int(delay.days + 0.5))
"""


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

folder_to_display = {v:k for k,v in displayname.items()}
folder_to_display.update({
    "United Kingdom": "UnitedKingdom",
    "Bosnia and Herzegovina": "BosniaandHerzegovina",
    "San Marino": "SanMarino",
            })

def display_to_folder(c):
    return folder_to_display.get(c, c)

fig, ax = plt.subplots(1, 1, figsize=(5, 5.3))


def Color(i):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    return colors[i % len(colors)]


def GlobColors():
    color = dict()
    for i, path in enumerate(
            sorted(glob(os.path.join(datafolder, "*", "intervals.json")))):
        folder = re.findall(os.path.join(datafolder, "(.*)", "intervals.json"),
                            path)[0]
        color[folder] = Color(i)
    return color

#color = {c: Color(i) for i, c in enumerate(sorted(csv['country']))}
color = GlobColors()

for i in range(len(csv)):
    date = pd.to_datetime(csv['date'][i])
    exact = pd.to_datetime(csv['exactDate'][i])

i = 0

last_day = pd.to_datetime(LAST_DAY)

yticks_y = []
yticks_label = []

myFmt = mdates.DateFormatter('%b %d')
ax.xaxis.set_major_formatter(myFmt)

argsort = np.argsort(csv['date'])

ax.get_yaxis().set_visible(False)

for i in range(len(csv)):
    ii = argsort[i]
    c = csv['country'][ii]
    date = pd.to_datetime(csv['date'][ii])
    exact = pd.to_datetime(csv['exactDate'][ii])

    y = -i

    yticks_y.append(y)
    yticks_label.append(c)

    cl = color[display_to_folder(c)]
    fcl = 'none' if pd.isna(exact) else cl
    p = ax.scatter(date, y, s=16, marker='o', facecolor=fcl, edgecolor=cl)

    if not pd.isna(exact):
        ax.scatter(exact, y, s=16, c=cl, marker='|')
        ax.plot([exact, date], [y, y], c=cl)
    ax.annotate(displayname.get(c, c),
                xy=(date, y),
                fontsize=7,
                xytext=(10, 0),
                textcoords='offset points',
                va='center')

ax.set_xticks(
    list(
        map(pd.to_datetime, [
            "2020-03-01",
            "2020-03-15",
            "2020-04-01",
            "2020-04-15",
            "2020-05-01",
            "2020-05-15",
            "2020-06-01",
        ])))

ax.set_yticks(yticks_y)
ax.set_yticklabels(yticks_label, fontsize=7)

fig.tight_layout()
fig.savefig("scatter_tact_actual.pdf")
