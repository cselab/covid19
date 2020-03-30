#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

from epidemics.tools.tools import import_from, load_file

import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--dataFolder', '-df', help='Load model data from this directory.')
args = parser.parse_args()






file = os.path.join(args.dataFolder,'initials.pickle')

x = load_file(file,'','pickle')

model = import_from( x['moduleName'], 'epModel' )

file = os.path.join(args.dataFolder,'state.pickle')

a = model( file )


a.plot_intervals()
