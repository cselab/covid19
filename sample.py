#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import argparse
import copy
from epidemics.tools.tools import import_from

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--compModel', '-cm', default='sir.nbin', help='The computational mode.')
parser.add_argument('--dataFolder', '-df', default='data/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--nSamples', '-ns', type=int, default=2000, help='Number of samples for TMCMC.')
parser.add_argument('--nGenerations', '-ng', type=int, default=20, help='Maximum number of generations for CMA-ES.')
parser.add_argument('--nThreads', '-nt', type=int, default=12, help='Number of threads.')
parser.add_argument('--nPropagation', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--futureDays', '-fd', type=int, default=2, help='Propagate that many days in future, after the time of observation of the last data.')
parser.add_argument('--nValidation', '-nv', type=int, default=0, help='Use that many data from the end of the data list to validate the prediction.')
parser.add_argument('--percentages', '-p', nargs='+', type=float, default=[0.5, 0.95, 0.99], help='Percentages for confidence intervals.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
parser.add_argument('--preprocess','-pre',type=bool,default=False,help='Preprocess infection data')
parser.add_argument('--backend','-b',default='numpy',help='numpy/scipy or torch')
parser.add_argument('--it_per_day','-it',default=10,help='number of iterations per day when using the pytorch ODE solver')

args = parser.parse_args()


x = copy.deepcopy(args)
del x.compModel
del x.nSamples
del x.nGenerations

model_class = import_from( 'epidemics.' + args.compModel, 'Model')

a = model_class( **vars(x) )

a.sample( args.nSamples )

a.propagate()

a.save()

a.plot_intervals(ns=20)
