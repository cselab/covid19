#!/usr/bin/env python

import matplotlib.pyplot as plt
import json
import os
from glob import glob
import re
import numpy as np
import adjustText
import sys

import countries
import countrydata

def printerr(m ,end='\n'):
    sys.stderr.write(str(m) + end)
    sys.stderr.flush()

# path to folder with output from `request_country.py`
datafolder = "."
df = countrydata.CollectCountryData(datafolder)
df = countrydata.AppendColor(df)
df = countrydata.AppendInferred(df, datafolder)

fig, axes = plt.subplots(1, 2, figsize=(9, 4))

for before, parname, ax in zip([True, False], ['R0int', 'R0'], axes):
    printerr(parname)
    ax.axvline(x=1, color='black', linestyle=':', zorder=-10)
    i = 0
    texts = []
    for i, row in enumerate(df.itertuples()):
        tint = (row.tint_mean - row.startday).days
        xy = [getattr(row, parname + "_mean"), tint]
        p = ax.scatter(*xy, s=16, c=row.color)
        texts.append(ax.annotate(countries.ABBREV2[row.folder], xy=xy, fontsize=6))
    printerr("adjusting text locations... ", end="")
    adjustText.adjust_text(texts, lim=10, on_basemap=True, ax=ax)
    printerr("done")
    ax.set_ylim(0, 70)
    ax.text(0.2 if not before else 0.2,
            1.03,
            r"$R_0$ {:} intervention".format("before" if before else "after"),
            transform=ax.transAxes,
            fontsize=15)
    ax.set_ylabel('days from first cases to intervention')
fig.tight_layout()
fig.savefig("scatter_tact.pdf")
