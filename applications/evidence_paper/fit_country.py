#!/usr/bin/env python3
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

sys.path.append('../../')
sys.path.append('../../build')


import argparse
import copy
from epidemics.utils.misc import import_from

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dataFolder', '-df', default='data/test/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--lastDay', '-ld', default='2020-06-15', help='Last day of data sequence in format %Y-%m-%d./')
parser.add_argument('--silentPlot', '-sp', action='store_true', help='Close plot window after plot.')
parser.add_argument('--preprocess', '-pre', type=bool, default=False, help='Preprocessing.')

args = parser.parse_args()

obs = []
obs.append('infections')
obs.append('deaths')

x = copy.deepcopy(args)
x.observations=obs

model_class = import_from( 'epidemics.country.reparam.sird_int.nbin', 'Model' )

a = model_class( **vars(x) )

a.fit_curve()
