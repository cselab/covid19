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

from countries import LAST_DAY


def days_to_delta(t):
    return np.timedelta64(int(t + 0.5), 'D')


folder_to_fullname = {
    "RussianFederation": "Russia",
    "UnitedKingdom": "United Kingdom",
    "BosniaandHerzegovina": "Bosnia and Herzegovina",
    "NorthMacedonia": "North Macedonia",
    "CzechRepublic": "Czech Republic",
    "VaticanCity": "Vatican City",
    "SanMarino": "San Marino",
}


def display_to_folder(c):
    return folder_to_display.get(c, c)


fig, ax = plt.subplots(1, 1, figsize=(5, 5.3))


def Color(i):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    return colors[i % len(colors)]


def CollectCountryData(datafolder):
    """
    Returns `pandas.DataFrame()` with data on countries found in `datafolder`.
        folder: subfolder in `datafolder` containing `intervals.json`
        startday: day of first data point
    """
    v_folder = []
    v_color = []
    v_startday = []
    v_fullname = []
    v_displayname = []
    for i, path in enumerate(
            sorted(glob(os.path.join(datafolder, "*", "intervals.json")))):
        folder = re.findall(os.path.join(datafolder, "(.*)", "intervals.json"),
                            path)[0]
        v_folder.append(folder)
        v_color.append(Color(i))
        with open(path) as f:
            js = json.loads(f.read())
        v_startday.append(
            np.datetime64(LAST_DAY) - days_to_delta(max(js['x-data'])))
        fullname = folder_to_fullname.get(folder, folder)
        v_fullname.append(fullname)
        v_displayname.append(fullname)

    countrydata = pd.DataFrame({
        "folder": v_folder,
        "color": v_color,
        "startday": v_startday,
        "displayname": v_displayname,
        "fullname": v_fullname,
    })
    return countrydata

def AppendOfficalLockdown(countrydata, path="official_lockdowns.csv"):
    """
    Appends country data from CollectCountryData() by dates of official lockdown.
    Columns: `[..., official_lockdown]`
    """
    csv = pd.read_csv(path)
    v_date = []
    for c in countrydata['fullname']:
        date = csv['date'].loc[csv['country'] == c]
        if not date.empty:
            v_date.append(pd.to_datetime(date.values.min()))
        else:
            v_date.append(pd.NaT)
    countrydata['official_lockdown'] = v_date

def AppendInferred(countrydata, datafolder):
    """
    Appends country data from CollectCountryData() by inferred intervenion time.
    Columns: `[..., official_lockdown]`
    """
    v_mean = []
    v_low = []
    v_high = []
    for folder,startday in zip(countrydata['folder'], countrydata['startday']):
        samples_path = os.path.join(datafolder, folder, 'sample_params.dat')
        assert os.path.isfile(samples_path)
        samples = np.genfromtxt(samples_path, names=True)
        tint = samples['tint']
        mean = np.mean(tint)
        low = np.quantile(tint, 0.10)
        high = np.quantile(tint, 0.90)
        v_mean.append(startday + days_to_delta(mean))
        v_low.append(startday + days_to_delta(low))
        v_high.append(startday + days_to_delta(high))
    countrydata['tint_mean'] = v_mean
    countrydata['tint_low'] = v_low
    countrydata['tint_high'] = v_high


datafolder = "."
countrydata = CollectCountryData(datafolder)
AppendOfficalLockdown(countrydata)
AppendInferred(countrydata, datafolder)
print(countrydata)


yticks_y = []
yticks_label = []

myFmt = mdates.DateFormatter('%b %d')
ax.xaxis.set_major_formatter(myFmt)
ax.get_yaxis().set_visible(False)

for i,row in enumerate(countrydata.sort_values(by='tint_mean').itertuples()):
    y = -i
    color = row.color
    official = row.official_lockdown
    ax.scatter(row.tint_mean,
               y,
               s=16,
               marker='o',
               facecolor='none' if pd.isnull(official) else color,
               edgecolor=color)
    ax.plot([row.tint_low, row.tint_high], [y, y],
            c="black",
            lw=3,
            alpha=0.25,
            zorder=-10)

    yticks_y.append(y)
    yticks_label.append(row.displayname)

    if not pd.isnull(official):
        ax.scatter(official, y, s=16, c=color, marker='|')
        ax.plot([official, row.tint_mean], [y, y], c=color)
    ax.annotate(row.displayname,
                xy=(row.tint_mean, y),
                fontsize=7,
                xytext=(15, 0),
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
        ])))

ax.set_yticks(yticks_y)
ax.set_yticklabels(yticks_label, fontsize=7)

fig.tight_layout()
fig.savefig("scatter_tint_actual.pdf")

