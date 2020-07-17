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

from epidemics.utils.misc import printlog, abort, moving_average
import libepidemics
import json
import numpy as np


def get_mean_parameters(dataFolder):
    from glob import glob
    gens = glob(os.path.join(dataFolder, "_korali_samples", "gen*.json"))
    lastgen = sorted(gens)[-1]
    with open(lastgen) as f:
        js = json.load(f)
    # FIXME: should be mean over all samples,
    #        'Mean Theta' is weighted by likelihood
    return js['Solver']['Mean Theta']


parser = argparse.ArgumentParser()
aa = parser.add_argument
aa('--dataFolder', default='data', help='Save all results in this folder')
aa('--data', nargs='+', type=float, help='Total infected.')
aa('--dataDays',
   nargs='+',
   type=float,
   help='Days at which `data` is defined, defaults to range(len(data)).')
aa('--populationSize', type=int, help='Total population.')
aa('--nSamples', type=int, default=2000, help='Number of samples for TMCMC.')
aa('--nThreads', type=int, default=1, help='Number of threads.')
aa('--nPropagate',
   type=int,
   default=100,
   help='Number of points to evaluate the solution in the propagation phase.')
aa('--nIntervals',
   type=int,
   default=100,
   help='Number of points for computing '
   'means and credible intervals.')
aa('--futureDays',
   '-fd',
   type=int,
   default=2,
   help=
   'Propagate that many days in future, after the time of observation of the last data.'
   )
aa('--percentages',
   nargs='+',
   type=float,
   default=[0.5, 0.9],
   help='Percentages for confidence intervals.')
aa('--duration',
   type=float,
   default=10,
   help='Duration of applying the  intervention.')
aa('--silent', action='store_true', help='No output on screen.')
aa('--moving_average',
   type=int,
   default=0,
   help='Half-width of moving average window applied to data.')
aa('--infer_duration',
   action='store_true',
   help='Infer the duration of intervention from the data.')
aa('--configure',
   action='store_true',
   help='Configure the model and save in dataFolder.')
aa('--sample', action='store_true', help='Sample model saved in dataFolder')
aa('--propagate',
   action='store_true',
   help='Propagate model saved in dataFolder')
aa('--intervals',
   action='store_true',
   help='Compute intervals for model saved in dataFolder')
args = parser.parse_args()

dataFolder = os.path.join(os.path.abspath('.'), args.dataFolder) + '/'
statefile = os.path.join(dataFolder, 'state.pickle')

from model import Model

os.makedirs(dataFolder, exist_ok=True)

if not any([args.configure, args.sample, args.propagate, args.intervals]):
    args.configure = True
    args.sample = True
    args.propagate = True
    args.intervals = True

if args.configure:
    params_to_infer = ['R0', 'tint', 'kint']
    if args.infer_duration:
        params_to_infer.append('dint')

    params_prior = dict()

    params_fixed = {
        'dint': args.duration,
    }

    data = args.data
    data = moving_average(data, args.moving_average)

    assert args.populationSize is not None,\
            "--populationSize is required"

    kwargs = {
        'dataTotalInfected': data,
        'dataFolder': dataFolder,
        'populationSize': args.populationSize,
        'dataDays': args.dataDays,
        'params_to_infer': params_to_infer,
        'params_prior': params_prior,
        'params_fixed': params_fixed,
    }

    model = Model(**kwargs)
    model.save(statefile)
else:
    model = Model.load(statefile)

if args.sample:
    model.sample(args.nSamples)
    model.save(statefile)
else:
    model = Model.load(statefile)

if args.propagate:
    model.propagate(min(model.nSamples, args.nPropagate), args.futureDays)
    model.save(statefile)
else:
    model = Model.load(statefile)

if args.intervals:
    from epidemics.epidemics import load_param_samples
    percentages = args.percentages
    vv = dict()
    vv['Daily Infected'] = model.compute_intervals("Daily Incidence",
                                                   args.nIntervals,
                                                   percentages)
    vv['Total Infected'] = model.compute_intervals("Daily Incidence",
                                                   args.nIntervals,
                                                   percentages,
                                                   cumsum=True)

    js = dict()
    js['x-data'] = list(model.data['Model']['x-data'])
    js['y-data'] = list(np.cumsum(model.data['Model']['y-data']))
    js['Population Size'] = model.populationSize
    js['nSamples'] = model.nSamples
    js['percentages'] = percentages

    samples = load_param_samples(dataFolder)
    mean_params = model.substitute_inferred(
        [np.mean(samples[k]) for k in model.params_to_infer])
    js['mean_params'] = mean_params
    js['intervals_params'] = {
        k: [(p, m, m) for p in percentages]
        for k, m in mean_params.items()
    }
    js['intervals_params'].update({
        k: [
            (
                p,  #
                np.quantile(samples[k], 0.5 - p / 2),  #
                np.quantile(samples[k], 0.5 + p / 2),  #
            ) for p in percentages
        ]
        for k in samples.dtype.names
    })

    for k, v in vv.items():
        t, mean, median, intervals = v
        js['x-axis'] = list(t)
        r = dict()
        r['Intervals'] = [{
            "Percentage": float(q[0]),
            "Low Interval": list(q[1]),
            "High Interval": list(q[2]),
        } for q in intervals]
        r['Mean'] = list(mean)
        r['Median'] = list(median)
        js[k] = r

    fn = os.path.join(dataFolder, 'intervals.json')
    printlog(f'Save intervals in: {fn}')
    with open(fn, 'w') as f:
        json.dump(js, f, indent=2, sort_keys=True)
