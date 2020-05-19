#!/usr/bin/env python

import geopandas as gpd
import pandas as pd
import geoplot
import matplotlib.pyplot as plt
import json
import os


import mapclassify

# XXX path to folder with output from `request_country.py`
datafolder = "data"

def get_folder(country):
    return os.path.join(datafolder, country.replace(' ', ''))
    #return country.replace(' ', '')


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

geocountries = map(get_geocountry, countries)

g_R0 = dict()
g_lambda = dict()
g_lambda2 = dict()
for c in countries:
    with open(os.path.join(get_folder(c), "intervals.json")) as f:
        js = json.loads(f.read())
    gc = get_geocountry(c)
    p = js['mean_params']
    g_R0[gc] = p['R0']
    gamma = p['gamma']
    beta = p['R0'] * gamma
    beta2 = beta * p['kbeta']
    g_lambda[gc] = beta - gamma
    g_lambda2[gc] = beta2 - gamma

df_R0 = pd.DataFrame(g_R0.items(), columns=['name', 'R0'])
df_lambda = pd.DataFrame(g_lambda.items(), columns=['name', 'lambda'])
df_lambda2 = pd.DataFrame(g_lambda2.items(), columns=['name', 'lambda2'])


world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world = world.to_crs(epsg=3395)
shapes = world[world['name'].isin(geocountries)]

shapes = shapes.merge(df_R0, on='name')
shapes = shapes.merge(df_lambda, on='name')
shapes = shapes.merge(df_lambda2, on='name')

fig, [ax,ax2] = plt.subplots(2)
ax.set_axis_off()
ax2.set_axis_off()

shapes.plot(column='lambda',
            ax=ax,
            edgecolor='black',
            lw=0.1,
            cmap=plt.get_cmap('coolwarm'),
            vmin=0,
            vmax=0.3)

shapes.plot(column='lambda2',
            ax=ax2,
            edgecolor='black',
            lw=0.1,
            cmap=plt.get_cmap('coolwarm'),
            vmin=-0.07,
            vmax=0.05)

plt.savefig("a.pdf")
