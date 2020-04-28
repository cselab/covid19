#!/usr/bin/env python3
# Author: Martin Boden
# Date:   21/04/2020
# Email:  mboden@ethz.ch

import sys
import korali

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
    e["File Output"]["Path"] = phase_3_path+'/_korali_samples_phase_3/'

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    k["Conduit"]["Type"] = "Concurrent"
    k["Conduit"]["Concurrent Jobs"] = 12
    k.run(e)

def propagation(model,canton,phase_3_path):
    model = 'sir.altone_nbin'
    n_samples = 1
    params = {'dataFolder': phase_3_path,
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

    a.load_parameters(phase_3_path+'/_korali_samples_phase_3/')
    a.sample()
    a.propagate()
    a.save()
    a.plot_intervals()

if __name__ == "__main__":  

    model = 'sir_altone_nbin'
    regions = CANTON_LIST_SHORT

    for canton in regions:
        print('Processing {}'.format(canton))

        phase_1_path = 'data/'+model+'/phase_1_results/'+canton+'/'+model+'/_korali_samples/'
        phase_2_path = 'data/'+model+'/phase_2_results/_korali_samples/'
        phase_3_path = 'data/'+model+'/phase_3_results/'+canton+'/'+model+'/'

        # sampling(phase_1_path,phase_2_path,phase_3_path)
        propagation(model,canton,phase_3_path)







