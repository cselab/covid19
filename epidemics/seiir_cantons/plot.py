#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import json
import itertools
from matplotlib import animation
import matplotlib.colors
import os
import collections

def hide_axis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)

code_to_name = {
    'ZH':'Zürich',
    'BE':'Bern',
    'LU':'Luzern',
    'UR':'Uri',
    'SZ':'Schwyz',
    'OW':'Obwalden',
    'NW':'Nidwalden',
    'GL':'Glarus',
    'ZG':'Zug',
    'FR':'Fribourg',
    'SO':'Solothurn',
    'BS':'Basel-Stadt',
    'BL':'Basel-Landschaft',
    'SH':'Schaffhausen',
    'AR':'Appenzell Ausserrhoden',
    'AI':'Appenzell Innerrhoden',
    'SG':'St. Gallen',
    'GR':'Graubünden',
    'AG':'Aargau',
    'TG':'Thurgau',
    'TI':'Ticino',
    'VD':'Vaud',
    'VS':'Valais',
    'NE':'Neuchâtel',
    'GE':'Genève',
    'JU':'Jura',
}

name_to_code = {}
for code,name in code_to_name.items():
    name_to_code[name] = code

codes = code_to_name.keys()


class Renderer:
    def __init__(self, code_to_timeseries):
        '''
        code_to_timeseries[code]: `numpy.ndarray`, (N)
            Dictionary with timeseries for every canton code.
        '''
        basedir = os.path.dirname(__file__)
        with open(os.path.join(basedir, "data/home_work_people.json")) as f:
            home_work_people = json.load(f)
        fname = os.path.join(basedir, "data/canton_shapes.npy")
        d = np.load(fname, allow_pickle=True).item()

        # Compute shape centers.
        centers = {}
        for name, ss in d.items():
            for i,s in enumerate(ss):
                x, y = s
                centers[name_to_code[name]] = [x.mean(), y.mean()]
                break

        colors = {}
        colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i,code in enumerate(codes):
            colors[code] = np.array(matplotlib.colors.to_rgb(
                    colorcycle[i % len(colorcycle)]))
        self.colors = colors

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
            code = name_to_code[name]
            for i,s in enumerate(ss):
                x, y = s
                line, = ax.plot(x, y, marker=None, c='black', lw=0.5)
                fill, = ax.fill(x, y, alpha=0.25, c='white')
                fills[code].append(fill)

        # Draw labels.
        texts = dict()
        self.texts = texts
        for code in codes:
            text = ax.text(*centers[code], code, ha='center', zorder=10)
            texts[code] = text
            ax.scatter(*centers[code], color=colors[code], zorder=5)

        # Draw connections.
        max_people = np.max([v for vv in home_work_people.values() for v in vv.values()])
        for c_home, c_work in itertools.product(codes, repeat=2):
            x0, y0 = centers[c_home]
            x1, y1 = centers[c_work]
            n = home_work_people[c_home][c_work]
            alpha = np.clip(n / max_people * 100, 0, 1)
            lw = np.clip(n / max_people * 5, 1, 4)
            ax.plot([x0, x1], [y0, y1], color=colors[c_work], alpha=alpha, lw=lw)

        self.code_to_timeseries = code_to_timeseries
        self.timeseries_length = len(next(iter(code_to_timeseries.values())))

    def init(self):
        return [v for vv in self.fills.values() for v in vv] + list(self.texts.values())

    def update(self, timeidx=-1, silent=False):
        if timeidx == -1:
            timeidx = self.timeseries_length - 1
        if not silent:
            print("{:}/{:}".format(timeidx + 1, self.timeseries_length))
        for code,timeseries in self.code_to_timeseries.items():
            u = timeseries[timeidx]
            color = self.colors[code] * np.clip(u, 0, 1)
            self.texts[code].set_text("{:}\n{:.2f}".format(code, u))
            for fill in self.fills[code]:
                fill.set_color(color * u)
        return [v for vv in self.fills.values() for v in vv] + list(self.texts.values())

    def save_movie(self, frames=100, filename="a.mp4"):
        ani = animation.FuncAnimation(self.fig, self.update,
                frames=np.linspace(0, self.timeseries_length - 1, frames, dtype=int),
                init_func=self.init, blit=True)
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=2000)
        ani.save('a.mp4', writer=writer)

    def save_image(self, timeidx=-1, filename="a.png"):
        self.init()
        self.update(timeidx, silent=True)
        self.fig.savefig(filename)

if __name__ == "__main__":
    rend = Renderer({"ZH":[0, 0.5, 0.3]})
    #rend.save_movie(frames=3)
    rend.save_image()
