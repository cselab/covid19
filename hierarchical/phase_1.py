#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import argparse
import copy
import sys
sys.path.append('../')
from epidemics.tools.tools import import_from
from epidemics.data.files.canton_population import CANTON_LIST

# Phase 1 for Hierarchical Bayesian Inferance: Individual level fit

model = 'sir.altone_nbin'
regions = CANTON_LIST
n_samples = 2000
params = {'dataFolder': './data',
          'preprocess':True,
          'nThreads': 12,
          'nPropagation': 100,
          'futureDays': 2,
          'nValidation': 0,
          'percentages': [0.5, 0.95, 0.99],
          'silent': False}

for region in regions:

    params['country'] = region
    model_class = import_from( 'epidemics.' + model, 'Model')

    a = model_class(**params)

    a.sample(n_samples)

    a.propagate()

    a.save()

    a.plot_intervals()
