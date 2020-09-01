#!/usr/bin/env python3

# Importing computational model
import sys
import os
import korali

sys.path.append('../../')
from epidemics.cantons.data.canton_population import CANTON_LIST, CANTON_LIST_SHORT
import argparse

'''
    Hierearchical setup for epidemiology models 
'''

def create_folder(name):
    if not os.path.exists(name):
        os.makedirs(name)

def get_regions(regions):
    ## Select regions
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
    elif args.regions == 'g9':
        regions = ['canada','china','france','germany','italy','japan','russia','switzerland','uk','us']
        folder_name = '/g9'

    else:
        regions = args.regions
        folder_name = str(regions)
    return regions, folder_name

def get_variables(model):

    model_parts = model.split('.')
    model_name = model_parts[-2]
    stat_model = model_parts[-1]

    possible_variables = {
                'R0':       {'name':'R0','cond_type':'Normal',      'Mean':(0.0,10.0),'Std':(0.0,5.0)},     # Average reproduction number
                'D':        {'name':'D','cond_type':'Normal',       'Mean':(0.0,25.0),'Std':(0.0,10.0)},    # Average recovery period
                'Z':        {'name':'Z','cond_type':'Normal',       'Mean':(0.0,30.0),'Std':(0.0,10.0)},    # Average incubation period
                'eps':      {'name':'eps','cond_type':'Normal',     'Mean':(0.0,0.25),'Std':(0.0,0.2)},     # Death rate
                'mu':       {'name':'mu','cond_type':'Normal',      'Mean':(0.0,30.0),'Std':(0.0,10.0)},    # Reduction factor for unreported infected
                'alpha':    {'name':'alpha','cond_type':'Normal',   'Mean':(0.0,30.0),'Std':(0.0,10.0)},    # Percentage reported
                'tact':     {'name':'tact','cond_type':'Normal',    'Mean':(0.0,30.0),'Std':(0.0,10.0)},    # Intervention timne
                'dtact':    {'name':'dtact','cond_type':'Normal',   'Mean':(0.0,50.0),'Std':(0.0,10.0)},    # Intervention lengths
                'kbeta':    {'name':'kbeta','cond_type':'Normal',   'Mean':(0.0,1.0),'Std':(0.0,0.1)},      # R0 reduction factor
                'r':        {'name':'r','cond_type':'Normal',       'Mean':(0.0,15.0),'Std':(0.0,5)},      # Nbin stat param
                'Sigma':    {'name':'Sigma','cond_type':'Normal',   'Mean':(0.0,5.0),'Std':(0.0,2.5)}       # Nbin stat param
            }

    # Variables for each parameter
    if model_name == 'sird_int':
        variables = ['R0','D','eps','tact','dtact','kbeta']
    elif model_name == 'sird':
        variables = ['R0','D','eps']

    # Stat model
    if stat_model == 'nbin':
        variables.append('r')
    elif stat_model == 'tnrm':
        variables.append('Sigma')

    return [possible_variables[var] for var in variables ]


def run_phase_2(phase_1_paths,phase_2_path,variables):

    # Problem
    e = korali.Experiment()
    e["Problem"]["Type"] = "Hierarchical/Psi"

    for i in range(len(phase_1_paths)):
        print(phase_1_paths[i])
        subProblem = korali.Experiment()
        subProblem.loadState(phase_1_paths[i])
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

    #Solver
    # e["Solver"]["Type"] = "Sampler/TMCMC"
    # e["Solver"]["Population Size"] = 2000
    # e["Solver"]["Default Burn In"] = 3;
    # e["Solver"]["Target Coefficient Of Variation"] = 0.6
    # e["Solver"]["Covariance Scaling"] = 0.01

    e["Solver"]["Type"] = "Sampler/Nested"
    e["Solver"]["Resampling Method"] = "Multi Ellipse"
    e["Solver"]["Number Live Points"] = 1500
    e["Solver"]["Proposal Update Frequency"] = 1500
    e["Solver"]["Ellipsoidal Scaling"] = 1.10
    batch = 12
    e["Solver"]["Batch Size"] = batch
 
    e["Solver"]["Termination Criteria"]["Max Generations"] = 1e9
    e["Solver"]["Termination Criteria"]["Min Log Evidence Delta"] = 0.1
    e["Solver"]["Termination Criteria"]["Max Effective Sample Size"] = 25000

    e["Console Output"]["Verbosity"] = "Detailed"
    e["File Output"]["Path"] = phase_2_path
    e["File Output"]["Frequency"] = 5000
    create_folder(phase_2_path)    

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    # k["Conduit"]["Type"] = "Concurrent"
    # k["Conduit"]["Concurrent Jobs"] = batch
    # print('Launching Korali')
    k.run(e)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='country.reparam.sir_int.tnrm', help='Model type')
    parser.add_argument('--phase_1_path', '-p', default='./data/country.reparam.sir_int.tnrm/phase_1_results', help='Model type')
    parser.add_argument('--regions', '-r', default='all', help='Model type')
    parser.add_argument('--output', '-o', default='same', help='output path')

    args = parser.parse_args()
    model = args.model

    ## Get regions and models
    regions,folder_name = get_regions(args.regions)
    variables = get_variables(args.model)

    if args.output == 'same':
        output_path = args.phase_1_path
    else:
        output_path = args.output


    ## Paths
    phase_1_data = [args.phase_1_path+'/'+region+'/'+model+'/_korali_samples/latest' for region in regions]
    phase_2_path = output_path + '/_hierarchical/'+model+'/'+folder_name+'/phase_2_results/_korali_samples'
    print(phase_2_path)

    run_phase_2(phase_1_data,phase_2_path,variables)
