#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import sys
import korali
import argparse

sys.path.append('../../')
from epidemics.tools.tools import import_from
from epidemics.data.files.canton_population import CANTON_LIST, CANTON_LIST_SHORT

def sampling(phase_1_path,phase_2_path,phase_3_path):

    e = korali.Experiment()

    e["Problem"]["Type"]  = "Hierarchical/Theta"
    print(phase_1_path)
    e["Problem"]["Theta Problem Path"] = phase_1_path
    e["Problem"]["Psi Problem Path"] = phase_2_path

    e["Solver"]["Type"] = "TMCMC"
    e["Solver"]["Population Size"] = 2000
    e["Solver"]["Termination Criteria"]["Max Generations"] = 30
    e["Solver"]["Default Burn In"] = 1;
    e["Solver"]["Target Coefficient Of Variation"] = 0.6

    e["Console Output"]["Verbosity"] = "Detailed"
    e["File Output"]["Path"] = phase_3_path+'/_korali_samples/'

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    k["Conduit"]["Type"] = "Concurrent"
    k["Conduit"]["Concurrent Jobs"] = 12
    k.run(e)

def propagation(plot_path,model,canton,phase_3_path):
    n_samples = 1

    dataFolder = plot_path+model+'/phase_3_results/'
    params = {'dataFolder': dataFolder,
              'preprocess':True,
              'nThreads': 12,
              'nPropagation': 100,
              'futureDays': 2,
              'nValidation': 0,
              'percentages': [0.5, 0.95, 0.99],
              'silent': False}

    model_class = import_from( 'epidemics.' + model, 'Model')
    params['country'] = canton
    a = model_class(**params)

    a.load_parameters(phase_3_path+'/_korali_samples/')
    a.propagate()
    a.save()
    a.plot_intervals(ns=20)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='sir_altone_nbin', help='Model type')
    parser.add_argument('--regions', '-r', default='cantons', help='Model type')
    parser.add_argument('--dir', '-dir', default='./data/', help='Model type')

    args = parser.parse_args()

    model = args.model

    if args.regions == 'cantons':
        regions = CANTON_LIST
    elif args.regions == 'cantons_short':
        regions = CANTON_LIST_SHORT

    for canton in regions:
        print('Processing {}'.format(canton))

        phase_1_path = args.dir+model+'/phase_1_results/'+canton+'/'+model+'/_korali_samples/'
        phase_2_path = args.dir+model+'/phase_2_results/_korali_samples/'
        phase_3_path = args.dir+model+'/phase_3_results/'+canton+'/'+model+'/'

        sampling(phase_1_path,phase_2_path,phase_3_path)
        propagation(args.dir,model,canton,phase_3_path)







