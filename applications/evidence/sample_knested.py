#!/usr/bin/env python3
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

sys.path.append('../../')
sys.path.append('../../build')


import argparse
import copy
from epidemics.tools.tools import import_from
sys.path.append('../../build')

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--compModel', '-cm', default='country.sir.tnrm', help='The computational model.')
parser.add_argument('--dataFolder', '-df', default='data/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--nSamples', '-ns', type=int, default=2000, help='Number of Live Samples.')
parser.add_argument('--nPropagation', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--nGenerations', '-ng', type=int, default=20, help='Maximum number of generations for CMA-ES.')
parser.add_argument('--futureDays', '-fd', type=int, default=2, help='Propagate that many days in future, after the time of observation of the last data.')
parser.add_argument('--nValidation', '-nv', type=int, default=0, help='Use that many data from the end of the data list to validate the prediction.')
parser.add_argument('--percentages', '-p', nargs='+', type=float, default=[0.5, 0.95, 0.99], help='Percentages for confidence intervals.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
parser.add_argument('--sampler', '-sa', default='TMCMC', help='Choose sampler TMCMC or mTMCMC')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--preprocess', '-pre', type=bool, default=False, help='Preprocessing.')
parser.add_argument('--up_to_int', '-utint', type=bool, default=False, help='Use only data before intervention')

args = parser.parse_args()


x = copy.deepcopy(args)
del x.compModel
del x.nSamples
del x.nPropagation
del x.nGenerations

model_class = import_from( 'epidemics.' + args.compModel, 'Model')

a = model_class( **vars(x) )

a.sample_knested(nLiveSamples=args.nSamples,freq=args.nSamples)

a.propagate( args.nPropagation )

a.plot_intervals()
