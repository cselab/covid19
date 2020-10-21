#!/usr/bin/env python3
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

sys.path.append('../../')
sys.path.append('../../build')


import argparse
import copy
from epidemics.utils.misc import import_from, put_comment
sys.path.append('../../build')

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--compModel', '-cm', default='country.reparam.sir_int.tnrm', help='The computational model.')
parser.add_argument('--dataFolder', '-df', default='data/test/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--nSamples', '-ns', type=int, default=1500, help='Number of Live Samples.')
parser.add_argument('--dLogz', '-dlz', type=float, default=0.1, help='Remaining Log Evidence threshold.')
parser.add_argument('--batchSize', '-bs', type=float, default=1, help='Number of samples evaluated at each iteration.')
parser.add_argument('--nPropagation', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--nGenerations', '-ng', type=int, default=20, help='Maximum number of generations.')
parser.add_argument('--futureDays', '-fd', type=int, default=2, help='Propagate that many days in future, after the time of observation of the last data.')
parser.add_argument('--nValidation', '-nv', type=int, default=0, help='Use that many data from the end of the data list to validate the prediction.')
parser.add_argument('--percentages', '-p', nargs='+', type=float, default=[0.5, 0.95, 0.99], help='Percentages for confidence intervals.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--preprocess', '-pre', type=bool, default=False, help='Preprocessing.')
parser.add_argument('--up_to_int', '-utint', type=bool, default=False, help='Use only data before intervention')
parser.add_argument('--useIntervention', '-uint', action='store_true', help='Add intervention start to tact')
parser.add_argument('--plotMeanMedian', dest='plotMeanMedian', action='store_true', default=False, help='Plot mean and median of states.')
parser.add_argument('--useInfections', '-ui', action='store_true', help='Use infections to fit data.')
parser.add_argument('--useInformedPriors', '-uip', action='store_true', help='Use informed priors on D, Z, Zp and Y.')
parser.add_argument('--useDeaths', '-ud', action='store_true', help='Use deaths to fit data.')
parser.add_argument('--test', action='store_true', help="Test run. Not everything is tested.")
parser.add_argument('--synthetic', '-syn', action='store_true', required=False, help='Run with sunthetic data (file must be provided).') 
parser.add_argument('--dataFile', '-dat', required=False, type=str, help='Datafile for synthetic data.') 

args = parser.parse_args()
obs = []
if args.useInfections:
    obs.append('infections')
if args.useDeaths:
    obs.append('deaths')

x = copy.deepcopy(args)
x.observations=obs
del x.compModel
del x.nSamples
del x.batchSize
del x.dLogz
del x.nPropagation
del x.nGenerations
del x.useInfections
del x.useDeaths
del x.test

model_class = import_from( 'epidemics.' + args.compModel, 'Model')

a = model_class( **vars(x) )

a.sample_knested(nLiveSamples=args.nSamples, freq=args.nSamples, dlogz=args.dLogz, 
                 batch=args.batchSize, maxiter=(5 if args.test else 1e9))

if not args.test:
    a.propagate( args.nPropagation )

    a.plot_intervals()
