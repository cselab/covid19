#!/usr/bin/env python3

import json
import numpy as np
import argparse
import os
from datetime import datetime

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


def linear_trans(u0, u1, t, tc, teps):
    """
    Linear transition from u0 to u1 in interval `tc - teps < t < tc + teps`.
    """
    t0 = tc - teps
    t1 = tc + teps
    return u0 if t <= t0 else u1 if t >= t1 else \
        u0 + (u1 - u0) * (t - t0) / (t1 - t0)


intervention_trans = linear_trans

import matplotlib.dates as mdates

parser = argparse.ArgumentParser()
parser.add_argument('dataDir', help="Path to 'intervals.json' and '_korali_samples'.")
parser.add_argument('--title', type=str, default="", help="Title for the figure")
parser.add_argument('--output_dir',default=".", help="Directory for output images.")
args = parser.parse_args()

dirFiles = os.listdir(args.dataDir+'_korali_samples/')

intervalFile = args.dataDir + '/intervals.json'
samplesFile = args.dataDir + '_korali_samples/' + dirFiles[-1]

with open(intervalFile) as f:
    d = json.load(f)

with open(samplesFile) as f:
    samples = json.load(f)


fig = plt.figure(figsize=(6, 9))
ax, ax2, ax3 = fig.subplots(3)

xdata = np.array(d['x-data']).astype(float)

t = np.array(d['x-axis']).astype(float)
params = d['mean_params']
tact = params['tact']
dtact = params['dtact']
gamma = params['gamma']
R0 = params['R0']
kbeta = params['kbeta']


def days_to_delta(t):
    return (t * 24 * 3600).astype('timedelta64[s]')


day0 = np.datetime64('2020-05-17') - days_to_delta(max(xdata))
days = day0 + days_to_delta(t)

daysdata = day0 + days_to_delta(xdata)

lambdaeff = np.zeros_like(t)
myFmt = mdates.DateFormatter('%b %d')
ax.xaxis.set_major_formatter(myFmt)
ax2.xaxis.set_major_formatter(myFmt)
ax3.xaxis.set_major_formatter(myFmt)


db    = np.array(samples['Results']['Sample Database']);
R0    = db[:,0]
tact  = db[:,1]
kbeta = db[:,2]

Nt = len(t)
Ns = R0.shape[0]

lambdaeff = np.zeros((Nt,Ns))
for i in range(Nt):
  for j in range(Ns):
    beta =  R0[j]* intervention_trans(1., kbeta[j], t[i], tact[j], dtact * 0.5) * gamma
    lambdaeff[i,j] = beta - gamma


mean = np.mean(lambdaeff,axis=1)
ax.plot(days, mean)

median = np.quantile(lambdaeff, 0.5, axis=1)
q1 = np.quantile(lambdaeff, 0.05, axis=1)
q2 = np.quantile(lambdaeff, 0.95, axis=1)
ax.fill_between(days, q1, q2, alpha=0.5)


ax.axhline(y=0, color='black', linestyle=':')
ax.set_ylabel(r'Effective growth rate $\lambda^\ast(t)$')

ydata = np.array(d['y-data']).astype(float)

mean = d['Daily Infected']['Mean']
median = d['Daily Infected']['Median']
low = d['Daily Infected']['Intervals'][0]['Low Interval']
high = d['Daily Infected']['Intervals'][0]['High Interval']
ax2.plot(days, median)
ax2.fill_between(days, low, high, alpha=0.5)
ax2.scatter(daysdata, np.diff(ydata, prepend=0), c='black', s=4)
ax2.set_ylabel('Daily new reported')
ax2.set_xlim(left=days.min())
ax2.set_ylim(bottom=0)

mean = d['Total Infected']['Mean']
median = d['Total Infected']['Median']
low = d['Total Infected']['Intervals'][0]['Low Interval']
high = d['Total Infected']['Intervals'][0]['High Interval']
ax3.plot(days, median)
ax3.fill_between(days, low, high, alpha=0.5)
ax3.scatter(daysdata, ydata, c='black', s=4)
ax3.set_ylabel('Total new reported')
ax3.set_xlim(left=days.min())
ax3.set_ylim(bottom=0)

if args.title:
    ax.set_title(args.title)

os.makedirs(args.output_dir, exist_ok=True)
p = os.path.join(args.output_dir, "total.pdf")
print(p)
fig.tight_layout()
fig.savefig(p)
