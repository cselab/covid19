#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import argparse
import copy
import sys
# sys.path.append('../')
from epidemics.tools.tools import import_from

# Phase 1 for Hierarchical Bayesian Inferance: Individual level fit

model = 'sir.altone_nbin'
data_folder = './data'

countries = [   'switzerland',
                'spain',
                'italy',
                'france',
                'germany',
                'netherlands',
                'uk',
            ]

# countries = [   'VD',
#                 'ZH',
#             ]

n_samples = 2000

for country in countries:

    params = {'dataFolder': data_folder,
              'country': country,
              'nThreads': 1,
              'nPropagation': 100,
              'futureDays': 2,
              'nValidation': 0,
              'percentages': [0.5, 0.95, 0.99],
              'silent': False}

    model_class = import_from( 'epidemics.' + model, 'Model')

    a = model_class(**params)

    a.sample(n_samples)

    a.propagate()

    a.save()

    a.plot_intervals()
