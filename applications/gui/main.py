#!/usr/bin/env python3

# Created by Petr Karnakov on 25.05.2020
# Copyright 2020 ETH Zurich

"""
Backend for GUI https://cse-lab.ethz.ch/coronavirus/

Implementation of `korali-apps:5.coronavirus/main.py`
using models in `epidemics`.
"""

import os
import argparse

import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'build'))

from epidemics.tools.tools import printlog, abort, moving_average
import libepidemics

parser = argparse.ArgumentParser()
parser.add_argument('--dataFolder', '-df', default='data', help='Save all results in this folder')
parser.add_argument('--data', nargs='+', type=float, help='Total infected.')
parser.add_argument('--dataDays', nargs='+', type=float, help='Days at which `data` is defined, defaults to `range(len(data))`.')
parser.add_argument('--populationSize', '-ps', type=int, default=80000, help='Total population.')
parser.add_argument('--nSamples', '-ns', type=int, default=2000, help='Number of samples for TMCMC.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--nPropagation', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--futureDays', '-fd', type=int, default=2, help='Propagate that many days in future, after the time of observation of the last data.')
parser.add_argument('--validateData', '-vd', type=int, default=0, help='Use that many data from the end of the data list to validate the prediction.')
parser.add_argument('--percentages', '-p', nargs='+', type=float, default=[0.5, 0.9], help='Percentages for confidence intervals.')
parser.add_argument('--duration', type=float, default=10, help='Duration of applying the  intervention.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--moving_average', type=int, default=0, help='Half-width of moving average window applied to data.')
parser.add_argument('--infer_duration', action='store_true', help='Infer the duration of intervention from the data.')
parser.add_argument('--configure', action='store_true', help='Configure the model and save in dataFolder.')
parser.add_argument('--sample', action='store_true', help='Sample model saved in dataFolder')
parser.add_argument('--propagate', action='store_true', help='Propagate model saved in dataFolder')
parser.add_argument('--intervals', action='store_true', help='Compute intervals for model saved in dataFolder')
args = parser.parse_args()

args.dataFolder = os.path.join(os.path.abspath('.'), args.dataFolder) + '/'
statefile = os.path.join(args.dataFolder, 'state.pickle')

from model import Model

if not any([args.configure, args.sample, args.propagate, args.intervals]):
    args.congiture = True
    args.sample = True
    args.propagate = True
    args.intervals = True

if args.configure:
    params_to_infer = ['R0', 'tint', 'kint']
    if args.infer_duration:
        params_to_infer.append('dint')

    params_fixed = {
        'dint': args.duration,
    }

    data = args.data
    data = data[:-args.validateData - 1]
    data = moving_average(data, args.moving_average)

    kwargs = {
        'dataTotalInfected': data,
        'populationSize': args.populationSize,
        'percentages': args.percentages,
        'dataDays': args.dataDays,
        'futureDays': args.futureDays,
        'params_to_infer': params_to_infer,
        'params_fixed': params_fixed,
        'nSamples': args.nSamples,
        'nPropagation': args.nPropagation,
    }

    model = Model(**kwargs)
    model.save(statefile)
else:
    model = Model.load(statefile)

if args.sample:
    model.sample()
    model.save(statefile)
else:
    model = Model.load(statefile)

if args.propagate:
    model.propagate()
    model.save(statefile)
else:
    model = Model.load(statefile)

# Compute credible intervals
#jskIntervals = compute_intervals(args, jskPropagation, jsData)
