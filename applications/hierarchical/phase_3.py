#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import sys
import korali
import argparse
import os
sys.path.append('../../')
sys.path.append('../../build')

from epidemics.utils.misc import import_from
from epidemics.data.files.canton_population import CANTON_LIST, CANTON_LIST_SHORT

def sampling(phase_1_path,phase_2_path,phase_3_path):

    e = korali.Experiment()
    theta = korali.Experiment()
    psi = korali.Experiment()

    theta.loadState(phase_1_path)
    psi.loadState(phase_2_path)


    e["Problem"]["Type"]  = "Hierarchical/Theta"
    e["Problem"]["Theta Experiment"] = theta
    e["Problem"]["Psi Experiment"] = psi

    e["Solver"]["Type"] = "Sampler/TMCMC"
    e["Solver"]["Population Size"] = 1000
    e["Solver"]["Default Burn In"] = 2
    e["Solver"]["Max Chain Length"] = 1
    e["Solver"]["Target Coefficient Of Variation"] = 0.6

    e["Console Output"]["Verbosity"] = "Detailed"
    e["File Output"]["Path"] = phase_3_path+'/_korali_samples/'

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    # k["Conduit"]["Type"] = "Concurrent"
    # k["Conduit"]["Concurrent Jobs"] = 12
    k.run(e)

def propagation(model,region,phase_3_path):
    n_samples = 1

    dataFolder = phase_3_path
    params = {'dataFolder': dataFolder,
              'preprocess':True,
              'nThreads': 12,
              'nPropagation': 100,
              'futureDays': 2,
              'nValidation': 0,
              'percentages': [0.5, 0.95, 0.99],
              'silent': False}

    model_class = import_from( 'epidemics.' + model, 'Model')
    params['country'] = region
    a = model_class(**params)

    a.load_parameters(phase_3_path+'/_korali_samples/')
    a.propagate()
    a.save()
    a.plot_intervals(ns=20)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='country.reparam.sir_int.tnrm', help='Model type')
    parser.add_argument('--regions', '-r', default='all', help='Model type')
    parser.add_argument('--phase_1_path', '-p', default='./data/country.reparam.sir_int.tnrm/phase_1_results', help='Model type')

    args = parser.parse_args()
    model = args.model

    if args.regions == 'all':
        regions = [region for region in os.listdir(args.phase_1_path) if not region.startswith('_')]
        folder_name = '/all_countries'
    elif args.regions == '/cantons':
        regions = CANTON_LIST
        folder_name = 'cantons'
    elif args.regions == '/cantons_short':
        regions = CANTON_LIST_SHORT
        folder_name = '/cantons_short'
    else:
        regions = args.regions
        folder_name = str(regions)

    phase_2_path = args.phase_1_path + '/_hierarchical/'+model+folder_name+'/phase_2_results/_korali_samples/latest'


    for region in regions:
        print('Processing {}'.format(region))

        phase_1_path = args.phase_1_path+'/'+region+'/'+model+'/_korali_samples/latest'
        phase_3_path = args.phase_1_path + '/_hierarchical/'+model+folder_name+'/phase_3_results/'+region

        print(phase_1_path,phase_2_path,phase_3_path)
        sampling(phase_1_path,phase_2_path,phase_3_path)

        # propagation(model,region,phase_3_path)







