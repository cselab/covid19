#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

from epidemics.tools.tools import import_from, load_file

import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--dataFolder', '-df', help='Load model data from this directory.', required=True)
args = parser.parse_args()


file = os.path.join(args.dataFolder,'initials.pickle')

x = load_file(file,'','pickle')

model_class = import_from( x['moduleName'], 'Model' )

file = os.path.join(args.dataFolder,'state.pickle')

a = model_class( file )

# a.propagate()
# a.compute_intervals()
# a.save()

a.plot_intervals()
