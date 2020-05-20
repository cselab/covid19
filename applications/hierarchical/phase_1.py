#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import sys
sys.path.append('../../')
from epidemics.tools.tools import import_from
from epidemics.data.files.canton_population import CANTON_LIST
import argparse
def run_phase_1(model,region,n_samples,params):

    model_class = import_from( 'epidemics.' + model, 'Model')
    params['country'] = region
    
    print('[Phase 1] Running {}'.format(region))
    a = model_class(**params)
    a.sample(n_samples)
    a.propagate()
    a.save()
    a.plot_intervals(ns=20)

if __name__ == "__main__": 

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='sir_int.nbin', help='Model type')
    parser.add_argument('--regions', '-r', default='cantons', help='Model type')
    parser.add_argument('--dir', '-dir', default='./data/', help='Model type')

    args = parser.parse_args()

    model = args.model

    if args.regions == 'cantons':
        regions = CANTON_LIST
    elif args.regions == 'cantons_short':
        regions = CANTON_LIST_SHORT

    n_samples = 2000
    params = {'dataFolder': args.dir'+model+'/phase_1_results/',
              'preprocess':True,
              'nThreads': 12,
              'nPropagation': 30,
              'futureDays': 2,
              'nValidation': 0,
              'percentages': [0.5, 0.95, 0.99],
              'silent': False}

    for region in regions:
        run_phase_1(model,region,n_samples,params)

