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
parser.add_argument('--lastDay', '-ld', default='2020-06-15', help='Last day of data sequence in format %Y-%m-%d./')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--preprocess', '-pre', type=bool, default=False, help='Preprocessing.')
parser.add_argument('--up_to_int', '-utint', type=bool, default=False, help='Use only data before intervention')
parser.add_argument('--useIntervention', '-uint', action='store_true', help='Use only data before intervention')
parser.add_argument('--useInfections', '-ui', action='store_true', help='Use infections to fit data.')
parser.add_argument('--useDeaths', '-ud', action='store_true', help='Use deaths to fit data.')
parser.add_argument('--test', action='store_true', help="Test run. Not everything is tested.")

args = parser.parse_args()
obs = []
if args.useInfections:
    obs.append('infections')
if args.useDeaths:
    obs.append('deaths')

x = copy.deepcopy(args)
x.observations=obs
del x.compModel
# add something for num samples
del x.useInfections
del x.useDeaths
del x.test


model_class = import_from( 'epidemics.' + args.compModel, 'Model')

a = model_class( **vars(x) )

a.sample_hmc()