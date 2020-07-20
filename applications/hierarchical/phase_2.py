#!/usr/bin/env python3

# Importing computational model
import sys
import os
import korali

sys.path.append('../../')
from epidemics.cantons.data.canton_population import CANTON_LIST, CANTON_LIST_SHORT
import argparse

'''
    Hierearchical setup for SIR models with negative binomial
'''


def run_phase_2_auto(phase_1_path,phase_2_path,variables):

    # Problem
    e = korali.Experiment()
    e["Problem"]["Type"] = "Hierarchical/Psi"
    for i in range(len(phase_1_path)):
        print(phase_1_path[i])
        subProblem = korali.Experiment()
        subProblem.loadState(phase_1_path[i])
        e["Problem"]["Sub Experiments"][i] = subProblem

    # Define conditionals
    e["Problem"]["Conditional Priors"] = ["Conditional "+str(var['name']) for var in variables ]

    for i, var in enumerate(variables):

        e["Distributions"][i]["Name"] = "Conditional "+str(var['name'])

        if var['cond_type'] == 'Normal':
            e["Distributions"][i]["Type"] = "Univariate/Normal"
            e["Distributions"][i]["Mean"] = var['name'] + ' Mean'
            e["Distributions"][i]["Standard Deviation"] = var['name'] + ' Std'
        else:
            print('not implemented yet')

    # Define hyperparameters
    distrib_counter = len(variables)
    i = 0
    for var in variables:
        cond_params = [ele for ele in list(var.keys()) if ele not in ['name','cond_type']] 

        for cond_param in cond_params:
            var_name = var['name'] + ' ' + cond_param
            e["Variables"][i]["Name"] = var_name
            e["Variables"][i]["Prior Distribution"] = 'Uniform ' + var_name

            j = distrib_counter+i  # offset to take into account the prior distributions
            e["Distributions"][j]["Name"] = 'Uniform '+ var_name
            e["Distributions"][j]["Type"] = "Univariate/Uniform"
            e["Distributions"][j]["Minimum"] = var[cond_param][0]
            e["Distributions"][j]["Maximum"] = var[cond_param][1]
            i += 1

    # Solver
    e["Solver"]["Type"] = "Sampler/TMCMC"
    e["Solver"]["Population Size"] = 2000
    e["Solver"]["Default Burn In"] = 3;
    e["Solver"]["Target Coefficient Of Variation"] = 0.6
    e["Solver"]["Covariance Scaling"] = 0.01

    e["Console Output"]["Verbosity"] = "Detailed"
    e["File Output"]["Path"] = phase_2_path

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    # k["Conduit"]["Type"] = "Concurrent"
    # k["Conduit"]["Concurrent Jobs"] = 12
    k.run(e)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='country.reparam.sir_int.tnrm', help='Model type')
    parser.add_argument('--phase_1_path', '-p', default='./data/country.reparam.sir_int.tnrm/phase_1_results', help='Model type')
    parser.add_argument('--regions', '-r', default='all', help='Model type')

    args = parser.parse_args()
    model = args.model

    if args.regions == 'all':
        regions = [region for region in os.listdir(args.phase_1_path) if not region.startswith('_')]
        print(regions)
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

    phase_1_data = [args.phase_1_path+'/'+region+'/'+model+'/_korali_samples/latest' for region in regions]

    phase_2_path = args.phase_1_path + '/_hierarchical/'+model+folder_name+'/phase_2_results/_korali_samples'


    variables = [   {'name':'R0','cond_type':'Normal',      'Mean':(0.1,10.0),'Std':(0.0,5.0)},
                    {'name':'D','cond_type':'Normal',       'Mean':(0.0,30.0),'Std':(0.0,10.0)},
                    {'name':'tact','cond_type':'Normal',    'Mean':(0.0,30.0),'Std':(0.0,10.0)},
                    {'name':'dtact','cond_type':'Normal',   'Mean':(0.0,50.0),'Std':(0.0,10.0)},
                    {'name':'kbeta','cond_type':'Normal',   'Mean':(0.0,1.0),'Std':(0.0,0.1)},
                    {'name':'[r]','cond_type':'Normal',     'Mean':(0.0,5.0),'Std':(0.0,2.5)}] 

    variables = [   {'name':'R0','cond_type':'Normal',      'Mean':(0.1,10.0),'Std':(0.0,5.0)},
                    {'name':'D','cond_type':'Normal',       'Mean':(0.0,30.0),'Std':(0.0,10.0)},
                    {'name':'[r]','cond_type':'Normal',     'Mean':(0.0,5.0),'Std':(0.0,2.5)}] 
    # run_phase_2(phase_1_data,phase_2_path)
    run_phase_2_auto(phase_1_data,phase_2_path,variables)
