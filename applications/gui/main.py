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

from epidemics.country.country import EpidemicsCountry
from epidemics.tools.tools import printlog, abort
import libepidemics

parser = argparse.ArgumentParser()
parser.add_argument('--dataFolder', '-df', default='data', help='Save all results in this folder')
parser.add_argument('--data', '-d', nargs='+', type=float, help='Infected population.')
parser.add_argument('--populationSize', '-ps', type=int, default=80000, help='Total population.')
parser.add_argument('--nSamples', '-ns', type=int, default=2000, help='Number of samples for TMCMC.')
parser.add_argument('--nSamplesPropagation', type=int, default=10, help='Number of samples for propagation to compute intervals.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--nPoints', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--futureDays', '-fd', type=int, default=2, help='Propagate that many days in future, after the time of observation of the last data.')
parser.add_argument('--validateData', '-vd', type=int, default=0, help='Use that many data from the end of the data list to validate the prediction.')
parser.add_argument('--percentages', '-p', nargs='+', type=float, default=[0.5, 0.9], help='Percentages for confidence intervals.')
parser.add_argument('--reduction', type=float, default=0.5, help='Reduction factor for R0 after intervention.')
parser.add_argument('--duration', type=float, default=10, help='Duration of applying the  intervention.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--moving_average', type=int, default=0, help='Half-width of moving average window applied to data.')
parser.add_argument('--infer_reduction', action='store_true', help='Infer the reduction factor for R0 after intervention from the data.')
parser.add_argument('--infer_duration', action='store_true', help='Infer the duration of intervention from the data.')
# parser.add_argument('--noSave', action='store_true', help='No intermediate files saved.')
args = parser.parse_args()
args.noSave = False  # until we fix the bug in korali

args.dataFolder = os.path.join(os.path.abspath('.'), args.dataFolder) + '/'

printlog(args)
abort(args)

# Download and save population data
#jsData = tools.download_data(args)

# Sample the parameters of the computational model
#jskSamples = sample_parameters(args, jsData)

# Propagate the uncertainty in the parameters
#jskPropagation = propagate_uncertainty(args, jskSamples, jsData)

# Compute credible intervals
#jskIntervals = compute_intervals(args, jskPropagation, jsData)
