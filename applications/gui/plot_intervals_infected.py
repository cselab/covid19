#!/usr/bin/env python3

# Created by Petr Karnakov on 29.05.2020
# Copyright 2020 ETH Zurich

import json
import numpy as np
import argparse
import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--dataFolder',
                    '-df',
                    default='data',
                    help='Load results from this folder')
parser.add_argument('--output_dir',
                    default=".",
                    help="Directory for output images.")
args = parser.parse_args()

with open(args.dataFolder + '/intervals.json') as f:
    d = json.load(f)

var = ['Total Infected', 'Daily Infected']

y_label = ['Total Infected', 'New Infected']

fig = plt.figure(figsize=(12, 8))
ax = fig.subplots(len(var))

ax[0].plot(d['x-data'],
           d['y-data'],
           'o',
           lw=2,
           label='infected (data)',
           color='black')
ax[0].set_ylim(0, np.max(d['y-data']) * 2)

ax[1].plot(d['x-data'],
           np.diff(d['y-data'], prepend=0),
           'o',
           lw=2,
           label='infected (data)',
           color='black')
ax[1].set_ylim(0, np.max(np.diff(d['y-data'])) * 2)

for k, y in enumerate(var):
    ax[k].plot(d['x-axis'],
               d[y]['Mean'],
               '-',
               lw=2,
               label='Mean',
               color='blue')
    ax[k].plot(d['x-axis'],
               d[y]['Median'],
               '-',
               lw=2,
               label='Median',
               color='black')

    d[y]['Intervals'].sort(key=lambda x: x['Percentage'], reverse=True)

    for x in d[y]['Intervals']:
        p1 = [max(k, 0) for k in x['Low Interval']]
        p2 = x['High Interval']
        p = 100 * x['Percentage']
        ax[k].fill_between(d['x-axis'],
                           p1,
                           p2,
                           alpha=0.5,
                           label=f' {p:.0f}% credible interval')

    ax[k].legend(loc='upper left')
    ax[k].set_ylabel(y_label[k])
    ax[k].grid()

ax[-1].set_xlabel('time in days')

os.makedirs(args.output_dir, exist_ok=True)
p = os.path.join(args.output_dir, "prediction.pdf")
print(p)
fig.tight_layout()
fig.savefig(p)
