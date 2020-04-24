#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from epidemics.tools.tools import load_model

from pathlib import Path
import argparse
import os

# required for load_model
from main import Model

dataFolder = Path("data")
f = dataFolder / 'cantons' / 'state.pickle'

assert os.path.isfile(f)

a = load_model(f)

# a.propagate()
# a.save()

a.plot_intervals()
a.plot_intervals(region=1)
