#!/usr/bin/env python3

from matplotlib import animation
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from pandas import Series

import collections
import itertools
import json
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data import DATA_CACHE_DIR, DATA_FILES_DIR
from epidemics.cantons.py.model import ModelData
from epidemics.data.swiss_cantons import NAME_TO_CODE, CODE_TO_NAME
from epidemics.tools.cache import cache, cache_to_file
import epidemics.data.swiss_municipalities as munic

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
            draw_Mij=True, draw_Cij=True):
        '''
        frame_callback: callable
            Function that takes Renderer and called before rendering a frame.
            It can use `set_values()` and `set_texts()` to update the state,
            and `get_frame()` and `get_max_frame()` to get current frame index.
        '''

        self.data = data
        self.frame_callback = frame_callback
        # FIXME: "code" vs "key" terminology, pick one.
        self.code_to_value = {}
        self.code_to_text = {}
        self.draw_zones = draw_zones

        fname = DATA_FILES_DIR / 'canton_shapes.npy'
        d = np.load(fname, allow_pickle=True).item()

        # Compute shape centers.
        centers = {}
        for name, ss in d.items():
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

        self.set_base_colors('red')

        dpi = 200
        fig, ax = plt.subplots(figsize=(1920 / dpi, 1080 / dpi), dpi=dpi)
        hide_axis(ax)
        ax.set_aspect('equal')
        fig.tight_layout()
        self.fig = fig
        self.ax = ax

        # Draw shapes.
        fills = collections.defaultdict(list)
        self.fills = fills
        for name, ss in d.items():
            code = NAME_TO_CODE[name]
            for i,s in enumerate(ss):
                x, y = s
                line, = ax.plot(x, y, marker=None, c='black', lw=0.25)
                fill, = ax.fill(x, y, alpha=0.25, c='white')
                fills[code].append(fill)

        # Draw zones
        if draw_zones:
            zone_to_canton, zone_to_geometry = munic.get_zones_info()
            zone_geoms = list(zone_to_geometry.values())
            gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(zone_geoms))
            gdf['names'] = list(zone_to_canton)
            gdf['cantons'] = list(zone_to_canton.values())
            gdf['values'] = np.zeros(gdf.shape[0])
            gdf_ax = gdf.plot(ax=self.ax, column='values', cmap='Greys', alpha=0.8)
            self.zone_gdf = gdf
            self.zone_gdf_ax = gdf_ax
            self.zone_to_canton = zone_to_canton

        # Draw labels.
        texts = dict()
        self.texts = texts
        for code in CODE_TO_NAME:
            xc, yc = centers[code]
            ax.text(xc, yc, code, ha='center', va='bottom', zorder=10,
                    color=[0,0,0])
            text = ax.text(
                    xc, yc - 1700,
                    '', ha='center', va='top', zorder=10, fontsize=7,
                    color=[0,0,0])
            texts[code] = text
            ax.scatter(xc, yc, color='black', s=8, zorder=5)

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
                alpha = np.clip(n / max_people * 20, 0, 0.5)
                lw = np.clip(n / max_people * 5, 0.5, 4)
                ax.plot([x0, x1], [y0, y1], color=color, alpha=alpha, lw=lw)

        if draw_Mij:
            _draw_connections(data.Mij, 'blue')
        if draw_Cij:
            _draw_connections(data.Cij, 'green')

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
            self.zone_gdf['values'] = zone_values

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
        return self.zone_gdf['values'].values

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

    def init(self):
        return [v for vv in self.fills.values() for v in vv] + list(self.texts.values())

    def update(self, frame=-1, silent=False):
        self.frame = frame
        self.frame_callback(self)

        if frame == -1:
            frame = self.max_frame
        if not silent:
            print("{:}/{:}".format(frame, self.max_frame))
        for code,value in self.code_to_value.items():
            color = self.base_colors[code]
            alpha = np.clip(value, 0, 1)
            for fill in self.fills[code]:
                fill.set_color(color)
                fill.set_alpha(alpha)
        for code in self.texts:
            if code in self.code_to_text:
                self.texts[code].set_text(str(self.code_to_text[code]))
        if self.draw_zones:
            self.zone_gdf.plot(ax=self.zone_gdf_ax, column='values', cmap='Greys', alpha=0.8)
        return [v for vv in self.fills.values() for v in vv] + list(self.texts.values())

    def save_movie(self, frames=100, filename="a.mp4", fps=15):
        self.max_frame = frames - 1
        ani = animation.FuncAnimation(self.fig, self.update,
                frames=list(range(frames)),
                init_func=self.init, blit=True)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=2000)
        ani.save(filename, writer=writer)

    def save_image(self, frame=-1, filename="a.png"):
        self.max_frame = 1
        self.init()
        self.update(self.max_frame, silent=True)
        self.fig.savefig(filename)

    def set_base_colors(self, code_to_rgb=None):
        if code_to_rgb is None:
            code_to_rgb = {}
            plt_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
            for i,code in enumerate(self.get_codes()):
                code_to_rgb[code] = plt_cycle[i % len(plt_cycle)]
        if isinstance(code_to_rgb, str):
            code_to_rgb = {code: code_to_rgb for code in self.get_codes()}
        self.base_colors = {}
        for code in code_to_rgb:
            self.base_colors[code] = np.array(matplotlib.colors.to_rgb(
                    code_to_rgb[code]))

if __name__ == "__main__":

    def frame_callback(rend):
        colors = dict()
        texts = dict()
        for i, c in enumerate(rend.get_codes()):
            colors[c] = np.sin(i + rend.get_frame() * 0.1) ** 2
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
            vv = np.random.rand(len(nn))
            rend.set_zone_values(vv)


    from epidemics.cantons.py.model import get_canton_model_data
    rend = Renderer(frame_callback, data=get_canton_model_data(),
            draw_zones=True)

    rend.save_image()
    rend.save_movie(frames=10, fps=5)
