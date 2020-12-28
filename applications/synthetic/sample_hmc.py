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
parser.add_argument('--compModel', '-cm', default='country.reparam.sir_int.nbin', help='The computational model.')
parser.add_argument('--dataFolder', '-df', default='data/test/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--nPropagation', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--preprocess', '-pre', type=bool, default=False, help='Preprocessing.')
parser.add_argument('--up_to_int', '-utint', type=bool, default=False, help='Use only data before intervention')
parser.add_argument('--useIntervention', '-uint', action='store_true', help='Use only data before intervention')
parser.add_argument('--useInfections', '-ui', action='store_true', help='Use infections to fit data.')
parser.add_argument('--useDeaths', '-ud', action='store_true', help='Use deaths to fit data.')
parser.add_argument('--test', action='store_true', help="Test run. Not everything is tested.")
parser.add_argument('--sampler', type=str, default='HMC')
parser.add_argument('--version', '-v', type=str, default='Euclidean')
parser.add_argument('--nSamples', '-ns', type=int, default=5000, help='Max iterations.')
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
del x.nPropagation
# add something for num samples
del x.useInfections
del x.useDeaths
del x.test
del x.nSamples


model_class = import_from( 'epidemics.' + args.compModel, 'Model')

a = model_class( **vars(x) )

a.sample_hmc( args.nSamples )
a.propagate( args.nPropagation )

a.plot_intervals()
