#!/usr/bin/env python

import geopandas as gpd
import pandas as pd
import geoplot
import matplotlib.pyplot as plt
import json
import os
from glob import glob
import re

import mapclassify

# XXX path to folder with output from `request_country.py`
datafolder = "."


def get_folder(country):
    return os.path.join(datafolder, country.replace(' ', ''))


def get_foldername(country):
    return country.replace(' ', '')


def get_geocountry(country):
    alias = {
        "Russian Federation": "Russia",
        "Czech Republic": "Czechia",
    }
    return alias.get(country, country)


countries = {
    "Russian Federation",
    "Ukraine",
    "France",
    "Spain",
    "Sweden",
    "Norway",
    "Germany",
    "Finland",
    "Poland",
    "Italy",
    "United Kingdom",
    "Romania",
    "Belarus",
    "Kazakhstan",
    "Greece",
    "Bulgaria",
    "Iceland",
    "Hungary",
    "Portugal",
    "Austria",
    "Czech Republic",
    "Serbia",
    "Ireland",
    "Lithuania",
    "Latvia",
    "Croatia",
    "Bosnia and Herzegovina",
    "Slovakia",
    "Estonia",
    "Denmark",
    "Switzerland",
    "Netherlands",
}


def geocountry_to_folder(gc):
    return gc.replace(' ', '')


geocountries = map(get_geocountry, countries)

world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world = world.to_crs(epsg=3395)

df_folder = world.name.copy()

for i, gc in world.name.iteritems():
    alias = {
        "S.Sudan": "SouthSudan",
        "DominicanRep.": "DominicanRepublic",
        "Russia": "RussianFederation",
        "CentralAfricanRep.": "CentralAfricanRepublic",
        "UnitedStatesofAmerica": "UnitedStates",
        "BosniaandHerz.": "BosniaandHerzegovina",
        "Czechia": "CzechRepublic",
    }
    folder = geocountry_to_folder(gc)
    df_folder[i] = alias.get(folder, folder)

world.insert(len(world.columns), 'folder', df_folder)

# folder name to parameter
f_R0 = dict()
f_lambda = dict()
f_lambda2 = dict()

folders = []
for path in glob(os.path.join(datafolder, "*", "intervals.json")):
    folder = re.findall(os.path.join(datafolder, "(.*)", "intervals.json"),
                        path)[0]
    with open(path) as f:
        js = json.loads(f.read())
    p = js['mean_params']
    gamma = p['gamma']
    beta = p['R0'] * gamma
    beta2 = beta * p['kbeta']
    folders.append(folder)
    f_R0[folder] = p['R0']
    f_lambda[folder] = beta - gamma
    f_lambda2[folder] = beta2 - gamma

shapes = world[world['folder'].isin(folders)]

#print("Missing folders:")
#print('\n'.join(f for f in folders if f not in shapes.folder.values))

df_R0 = pd.DataFrame(f_R0.items(), columns=['folder', 'R0'])
df_lambda = pd.DataFrame(f_lambda.items(), columns=['folder', 'lambda'])
df_lambda2 = pd.DataFrame(f_lambda2.items(), columns=['folder', 'lambda2'])

shapes = shapes.merge(df_R0, on='folder')
shapes = shapes.merge(df_lambda, on='folder')
shapes = shapes.merge(df_lambda2, on='folder')

#print("Found folders:")
#print('\n'.join(shapes.folder.values))

fig, [ax, ax2] = plt.subplots(2)
ax.set_axis_off()
ax2.set_axis_off()

ax.set_title(r"$\lambda*=\beta - \gamma$ before intervention")
shapes.plot(column='lambda',
            ax=ax,
            edgecolor='black',
            lw=0.1,
            cmap=plt.get_cmap('coolwarm'),
            vmin=0,
            vmax=0.3,
            legend=True)

ax2.set_title(r"$\lambda*=\beta - \gamma$ after intervention")
shapes.plot(column='lambda2',
            ax=ax2,
            edgecolor='black',
            lw=0.1,
            cmap=plt.get_cmap('coolwarm'),
            vmin=-0.1,
            vmax=0.05,
            legend=True)

plt.tight_layout()
plt.savefig("map.pdf")
