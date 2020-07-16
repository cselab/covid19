#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import argparse
import copy
from epidemics.utils.misc import import_from

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--compModel', '-cm', default='sir.altone_nbin', help='The computational mode.')
parser.add_argument('--dataFolder', '-df', default='data/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--lastDay', '-ld', default='2020-06-13', help='Last day of data sequence in format %Y-%m-%d./')
parser.add_argument('--preprocess', '-pre', type=bool, default=False, help='Preprocess data')
parser.add_argument('--nSamples', '-ns', type=int, default=16, help='Number of samples for CMAES.')
parser.add_argument('--nGenerations', '-ng', type=int, default=1000, help='Maximum number of generations for CMA-ES.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
parser.add_argument('--useInfections', '-ui', action='store_true', help='Use infections to fit data.')
parser.add_argument('--useDeaths', '-ud', action='store_true', help='Use deaths to fit data.')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
args = parser.parse_args()

x = copy.deepcopy(args)

del x.compModel
del x.nSamples
del x.nGenerations
del x.useInfections
del x.useDeaths

obs = []
if args.useInfections:
    obs.append('infections')
if args.useDeaths:
    obs.append('deaths')

x.observations=obs

model_class = import_from( 'epidemics.' + args.compModel, 'Model')

a = model_class( **vars(x) )

a.optimize( args.nSamples, args.nGenerations )

a.save()
