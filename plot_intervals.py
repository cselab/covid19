#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

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
# a.compute_intervals()
# a.save()

a.plot_intervals()
