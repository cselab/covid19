#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import sys
sys.path.append('../../')
from epidemics.tools.tools import import_from
from epidemics.data.files.canton_population import CANTON_LIST

def run_phase_1(model,region,n_samples,params):

    model_class = import_from( 'epidemics.' + model, 'Model')
    params['country'] = region
    
    print('[Phase 1] Running {}'.format(region))
    a = model_class(**params)
    a.sample(n_samples)
    a.propagate()
    a.save()
    a.plot_intervals()

if __name__ == "__main__": 

    CANTON_LIST = [ 'ZH',
                    'BE',
                    'LU',
                    'UR',
                    'SZ',
                    'OW',
                    'NW',
                    'GL',
                    'ZG',
                    'FR',
                    'SO',
                    'BS',
                    'BL',
                    'SH',
                    'AR',
                    'AI',
                    'SG',
                    'GR',
                    'AG',
                    'TG',
                    'TI',
                    'VD',
                    'VS',
                    'NE',
                    'GE',
                    'JU',
                ]

    model = 'sir_int.nbin'
    regions = CANTON_LIST
    n_samples = 2000
    params = {'dataFolder': './data/'+model.replace('.','_')+'/phase_1_results/',
              'preprocess':True,
              'nThreads': 12,
              'nPropagation': 100,
              'futureDays': 2,
              'nValidation': 0,
              'percentages': [0.5, 0.95, 0.99],
              'silent': False}

    for region in regions:
        run_phase_1(model,region,n_samples,params)

