import pandas as pd
import os
import re
import json
import numpy as np
from glob import glob

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


def CollectCountryData(datafolder):
    """
    Returns `pandas.DataFrame()` with data on countries found in `datafolder`.
        folder: subfolder in `datafolder` containing `intervals.json`
        startday: day of first data point
    """
    v_folder = []
    v_startday = []
    v_fullname = []
    v_displayname = []
    for i, path in enumerate(
            sorted(glob(os.path.join(datafolder, "*", "intervals.json")))):
        folder = re.findall(os.path.join(datafolder, "(.*)", "intervals.json"),
                            path)[0]
        v_folder.append(folder)
        with open(path) as f:
            js = json.loads(f.read())
        v_startday.append(
            np.datetime64(LAST_DAY) - days_to_delta(max(js['x-data'])))
        fullname = folder_to_fullname.get(folder, folder)
        v_fullname.append(fullname)
        v_displayname.append(fullname)

    df = pd.DataFrame({
        "folder": v_folder,
        "startday": v_startday,
        "displayname": v_displayname,
        "fullname": v_fullname,
    })
    return df


def AppendColor(df):
    """
    Appends country data from CollectCountryData() by color cycle from matplotlib.
    Columns: `[..., color]`
    """
    import matplotlib.pyplot as plt

    def Color(i):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        return colors[i % len(colors)]

    df['color'] = [Color(i) for i in range(len(df))]
    return df


def AppendOfficalLockdown(df, path="official_lockdowns.csv"):
    """
    Appends country data from CollectCountryData() by dates of official lockdown.
    Columns: `[..., official_lockdown]`
    """
    csv = pd.read_csv(path)
    v_date = []
    for c in df['fullname']:
        date = csv['date'].loc[csv['country'] == c]
        if not date.empty:
            v_date.append(pd.to_datetime(date.values.min()))
        else:
            v_date.append(pd.NaT)
    df['official_lockdown'] = v_date


def AppendInferred(df, datafolder):
    """
    Appends country data from CollectCountryData() by inferred intervenion time.
    Columns: `[..., official_lockdown]`
    """
    v_mean = []
    v_low = []
    v_high = []
    for folder, startday in zip(df['folder'], df['startday']):
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
    df['tint_mean'] = v_mean
    df['tint_low'] = v_low
    df['tint_high'] = v_high
