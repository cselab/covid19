#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from epidemics.tools.tools import load_model

from pathlib import Path
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--dataFolder', '-df', type=Path, help='Load model data from this directory.', required=True)
args = parser.parse_args()


file = args.dataFolder / 'state.pickle'

a = load_model( file )

# a.propagate()
# a.save()

a.plot_intervals()
