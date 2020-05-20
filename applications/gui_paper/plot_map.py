#!/usr/bin/env python

import geopandas as gpd
import pandas as pd
import geoplot
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
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


def geocountry_to_folder(gc):
    return gc.replace(' ', '')


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
        "Macedonia": "NorthMacedonia",
    }
    folder = geocountry_to_folder(gc)
    df_folder[i] = alias.get(folder, folder)

world.insert(len(world.columns), 'folder', df_folder)

# folder name to parameter
f_R0 = dict()
f_R0_int = dict()

folders = []
for path in glob(os.path.join(datafolder, "*", "intervals.json")):
    folder = re.findall(os.path.join(datafolder, "(.*)", "intervals.json"),
                        path)[0]
    with open(path) as f:
        js = json.loads(f.read())
    p = js['mean_params']
    gamma = p['gamma']
    kbeta = p['kbeta']
    folders.append(folder)
    f_R0[folder] = p['R0']
    f_R0_int[folder] = p['R0'] * kbeta

shapes = world[world['folder'].isin(folders)]

#print("Folders missing on the map:")
#print('\n'.join(f for f in folders if f not in shapes.folder.values))

df_R0 = pd.DataFrame(f_R0.items(), columns=['folder', 'R0'])
df_R0_int = pd.DataFrame(f_R0_int.items(), columns=['folder', 'R0_int'])

shapes = shapes.merge(df_R0, on='folder')
shapes = shapes.merge(df_R0_int, on='folder')

#print("Countries on the map:")
#print('\n'.join(world.name.values))

#print("Matched folders:")
#print('\n'.join(shapes.folder.values))

fig, [ax, ax2] = plt.subplots(1,2,figsize=(9,5.5))
#ax.set_axis_off()
#ax2.set_axis_off()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

switzerland = shapes[shapes.name == "Switzerland"]
bb = switzerland.total_bounds
ext = [5e6, 5e6]
shift = [1.4e6, 2e6]
xlim = [bb[0] - ext[0] + shift[0], bb[2] + ext[0] + shift[0]]
ylim = [bb[1] - ext[1] + shift[1], bb[3] + ext[1] + shift[1]]

class Norm(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

norm = Norm(vmin=0.5, vmax=2.5, vcenter=1)
cmap = plt.get_cmap('coolwarm')

ax.set_title(r"$R_0$ after intervention", fontsize=15)
shapes.plot(column='R0_int',
            ax=ax,
            edgecolor='black',
            lw=0.1,
            norm=norm,
            cmap=cmap,
            vmin=0,
            vmax=3)

ax.set_xlim(*xlim)
ax.set_ylim(*ylim)

ax2.set_title(r"$R_0$ before intervention", fontsize=15)
shapes.plot(column='R0',
            ax=ax2,
            edgecolor='black',
            lw=0.1,
            norm=norm,
            cmap=cmap,
            vmin=0,
            vmax=3)

ax2.set_xlim(*xlim)
ax2.set_ylim(*ylim)


fig.subplots_adjust(top=1,bottom=0.05, left=0.025, right=0.975, wspace=0.05)
cbar_ax = fig.add_axes([0.05, 0.07, 0.9, 0.05])
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
fig.colorbar(sm, cax=cbar_ax, cmap=cmap,norm=norm,orientation='horizontal',
        ticks=[0.5, 1, 1.5, 2, 2.5])

#fig.tight_layout()
fig.savefig("map.pdf")
