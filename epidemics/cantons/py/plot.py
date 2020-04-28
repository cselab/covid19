#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")

from matplotlib import animation
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from pandas import Series
import shapely

import collections
import itertools
import json
import os
import sys
import time

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data import DATA_CACHE_DIR, DATA_FILES_DIR
from epidemics.cantons.py.model import ModelData
from epidemics.data.swiss_cantons import NAME_TO_CODE, CODE_TO_NAME
from epidemics.tools.cache import cache, cache_to_file
import epidemics.data.swiss_municipalities as munic

def Log(msg):
    sys.stderr.write(str(msg) + "\n")


@cache
@cache_to_file(DATA_CACHE_DIR / 'zenodo_2017_centers.pickle')
def get_zone_centers():
    centers = {}
    zone_to_canton, zone_to_geometry = munic.get_zones_info()
    name_to_center = {}
    for name,geom in zone_to_geometry.items():
        name_to_center[name] = geom.centroid.coords[0]
    df = munic.get_name_and_population()
    key_to_name = Series(df.name.values,index=df.key).to_dict()
    for key,name in key_to_name.items():
        if name in name_to_center:
            centers[key] = list(name_to_center[name])
    return centers

def hide_axis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis='both', which='both', length=0)

code_to_center_shift = {
    'BE':(0,0),
    'LU':(0,0),
    'UR':(0,0),
    'SZ':(0,0),
    'OW':(0,-3),
    'NW':(2,2),
    'GL':(0,0),
    'ZG':(-2,-2),
    'FR':(5,-3),
    'SO':(-5,0),
    'BS':(-1,0.5),
    'BL':(9,-4),
    'SH':(-4,1),
    'AR':(-9,-6),
    'AI':(3,-3),
    'SG':(-7,-3),
    'GR':(0,0),
    'AG':(4,5),
    'TG':(5,7),
    'TI':(0,5),
    'VD':(-20,5),
    'VS':(0,0),
    'NE':(1,-3),
    'GE':(0,-2),
    'JU':(0,3),
        }

for key in code_to_center_shift:
    shift = code_to_center_shift[key]
    k = 1e3
    code_to_center_shift[key] = [shift[0] * k, shift[1] * k]


class Renderer:
    def __init__(self, frame_callback, data: ModelData, draw_zones=False,
            draw_Mij=True, draw_Cij=True, resolution=(1920,1080),
            airports=None):
        '''
        frame_callback: callable
            Function that takes Renderer.
            The function is called before rendering a frame,
            and `get_frame()` takes values between 0 and `get_max_frame()`.
            It can use `set_values()` and `set_texts()` to update the state,.
        '''

        self.data = data
        self.frame_callback = frame_callback
        # FIXME: "code" vs "key" terminology, pick one.
        self.code_to_value = {}
        self.code_to_text = {}
        self.draw_zones = draw_zones
        self.draw_Mij = draw_Mij
        self.draw_Cij = draw_Cij
        self.airports = airports

        fname = DATA_FILES_DIR / 'canton_shapes.npy'
        self.canton_shapes = np.load(fname, allow_pickle=True).item()

        # Compute shape centers.
        centers = {}
        for name, ss in self.canton_shapes.items():
            for i,s in enumerate(ss):
                x, y = s
                code = NAME_TO_CODE[name]
                centers[code] = [x.mean(), y.mean()]
                if code in code_to_center_shift:
                    shift = code_to_center_shift[code]
                    centers[code][0] += shift[0]
                    centers[code][1] += shift[1]
                break

        if draw_zones:
            centers.update(get_zone_centers())
        self.centers = centers

        resolution = np.array(resolution).astype(float)
        self.resolution = resolution
        dpi = resolution.min() / 5.
        figsize = resolution / dpi
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        self.fig = fig
        self.ax = ax
        self.last_frame_time = time.time()

    def set_values(self, code_to_value):
        '''
        code_to_value: `dict`
          Mapping from canton code to float between 0 and 1.
        '''
        self.code_to_value = code_to_value

    def set_zone_values(self, zone_values):
        '''
        code_to_value: `numpy.ndarray`, shape=self.zone_gdf.shape[0]
          Array of values between 0 and 1 to color the zones.
        '''
        if self.draw_zones:
            self.zone_values = zone_values

    def set_texts(self, code_to_text):
        '''
        code_to_text: `dict`
          Mapping from canton code to label text.
        '''
        self.code_to_text = code_to_text

    def get_values(self):
        return self.code_to_value

    def get_texts(self):
        return self.code_to_text

    def get_zone_names(self):
        return self.zone_gdf['names'].values

    def get_zone_values(self):
        return self.zone_values

    def get_zone_to_canton(self):
        return self.zone_to_canton

    def get_codes(self):
        return CODE_TO_NAME.keys()

    def get_frame(self):
        '''
        Returns current frame index between 0 and get_max_frame().
        '''
        return self.frame

    def get_max_frame(self):
        '''
        Returns maximum frame index.
        Set to `frames - 1` by `save_movie(frames)`.
        Set to 1 by `save_image()`.
        '''
        return self.max_frame

    def init_plot(self):
        '''
        Updates data and returns a list of artist to animate
        (would make an effect in case blit=True).
        '''
        self.ax.clear()
        ax = self.ax
        hide_axis(ax)
        ax.set_aspect('equal')
        self.fig.tight_layout()

        # Draw cantons.
        fills = collections.defaultdict(list)
        for name, ss in self.canton_shapes.items():
            code = NAME_TO_CODE[name]
            for i,s in enumerate(ss):
                x, y = s
                line, = ax.plot(x, y, marker=None, c='black', lw=0.25)
                fill, = ax.fill(x, y, alpha=0, c='red', lw=0)
                fills[code].append(fill)
        self.fills = fills

        # Draw zones (municipalities).
        if self.draw_zones:
            zone_to_canton, zone_to_geometry = munic.get_zones_info()
            zone_geoms = list(zone_to_geometry.values())
            gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(zone_geoms))
            gdf['names'] = list(zone_to_canton)
            gdf['cantons'] = list(zone_to_canton.values())
            self.zone_values = np.zeros(gdf.shape[0])
            self.zone_gdf = gdf
            self.zone_to_canton = zone_to_canton

            zone_fills = [[] for i in range(gdf.shape[0])]
            for i,geomvalue in enumerate(zip(gdf['geometry'], self.zone_values)):
                geom, value = geomvalue
                if isinstance(geom, shapely.geometry.MultiPolygon):
                    geom = list(geom)
                if isinstance(geom, shapely.geometry.Polygon):
                    geom = [geom]
                for g in geom:
                    poly = np.array(g.exterior.coords).T
                    fill, = ax.fill(*poly, alpha=0, c='black', lw=0, zorder=5)
                    zone_fills[i].append(fill)
            self.zone_fills = zone_fills


        data = self.data
        centers = self.centers
        def _draw_connections(matrix, color):
            max_people = np.max(matrix)
            if max_people == 0:
                return
            for c_home, c_work in itertools.product(data.region_keys, repeat=2):
                if c_home == c_work:
                    continue
                if c_home not in centers or c_work not in centers:
                    continue
                x0, y0 = centers[c_home]
                x1, y1 = centers[c_work]
                n = matrix[data.key_to_index[c_home], data.key_to_index[c_work]]
                alpha_min = 0.01
                alpha = np.clip(n / max_people * 20, alpha_min, 0.5)
                if alpha == alpha_min:
                    continue
                lw = np.clip(n / max_people * 5, 0.5, 4)
                ax.plot([x0, x1], [y0, y1], color=color, alpha=alpha, lw=lw)

        if self.draw_Mij:
            _draw_connections(data.Mij, 'blue')
        if self.draw_Cij:
            _draw_connections(data.Cij, 'green')

        # Draw labels.
        texts = dict()
        self.texts = texts
        for code in CODE_TO_NAME:
            xc, yc = self.centers[code]
            ax.text(xc, yc, code, ha='center', va='bottom', zorder=10,
                    color=[0,0,0])
            text = ax.text(
                    xc, yc - 1700,
                    '', ha='center', va='top', zorder=10, fontsize=7,
                    color=[0,0,0])
            texts[code] = text
            ax.scatter(xc, yc, color='black', s=8, zorder=5)

        if self.airports is not None:
            for a in self.airports:
                if a in centers:
                    marker="$\u2708$"
                    ax.scatter(*centers[a], marker=marker, zorder=6,
                            s=200, alpha=0.75, facecolor='green',
                            lw=0)

        return []

    def update_plot(self, frame=-1, silent=False):
        '''
        Updates data and returns a list of artist to update
        (would make an effect in case blit=True).
        '''
        plt.draw()
        self.frame = frame
        self.frame_callback(self)

        if frame == -1:
            frame = self.max_frame
        if not silent:
            time1 = time.time()
            dtime = time1 - self.last_frame_time
            Log("{:}/{:} {:.0f} ms".format(frame, self.max_frame, dtime * 1e3))
            self.last_frame_time = time1
        for code,value in self.code_to_value.items():
            alpha = np.clip(value, 0, 1) * 0.75
            for fill in self.fills[code]:
                fill.set_alpha(alpha)
        for code in self.texts:
            if code in self.code_to_text:
                self.texts[code].set_text(str(self.code_to_text[code]))
        if self.draw_zones:
            for fills,value in zip(self.zone_fills, self.zone_values):
                alpha = np.clip(value, 0, 1) * 0.25
                for fill in fills:
                    fill.set_alpha(alpha)
        return []

    def save_movie(self, frames=100, filename="a.mp4", fps=15,):
        self.max_frame = frames - 1
        ani = animation.FuncAnimation(self.fig, self.update_plot,
                frames=list(range(frames)),
                init_func=self.init_plot, blit=False)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=2000)
        ani.save(filename, writer=writer)

    def run_interactive(self, frames=100, fps=15):
        self.max_frame = frames - 1
        ani = animation.FuncAnimation(self.fig, self.update_plot,
                frames=list(range(frames)),
                init_func=self.init_plot, blit=False, interval=1000. / fps)
        plt.show()

    def save_image(self, frame=-1, frames=1, filename="a.png"):
        self.max_frame = frames - 1
        self.init_plot()
        self.update_plot(self.max_frame, silent=True)
        self.fig.savefig(filename)

def example():
    def frame_callback(rend):
        colors = dict()
        texts = dict()
        for i, c in enumerate(rend.get_codes()):
            colors[c] = np.sin(i + rend.get_frame() * 5 / rend.get_max_frame()) ** 2
            texts[c] = "{:},{:}".format(rend.get_frame(), i)
        rend.set_values(colors)
        rend.set_texts(texts)
        if rend.draw_zones:
            nn = rend.get_zone_names()
            sel_nn = nn[::100]
            vv = rend.get_zone_values()
            ii = []
            for i,n in enumerate(nn):
                if n in sel_nn:
                    ii.append(i)
            for i in ii:
                c = rend.get_zone_to_canton()[nn[i]]
                if c in colors:
                    vv[i] = colors[c]
            vv = np.cos(np.arange(len(nn)) + rend.get_frame() * 10 / rend.get_max_frame()) ** 2
            rend.set_zone_values(vv)


    from epidemics.cantons.py.model import get_canton_model_data
    rend = Renderer(frame_callback, data=get_canton_model_data(),
            draw_zones=True, resolution=(720,480))

    #rend.run_interactive(frames=50)
    rend.save_image()
    rend.save_movie(frames=10, fps=10)

if __name__ == "__main__":
    example()
