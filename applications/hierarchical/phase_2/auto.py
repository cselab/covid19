#!/usr/bin/env python3

# Importing computational model
import sys
import os
import korali

sys.path.append('../../')
from epidemics.data.files.canton_population import CANTON_LIST, CANTON_LIST_SHORT
import argparse

'''
    Hierearchical setup for SIR models with negative binomial
'''

    variables = [   {'name':'R0','cond_type':'Normal','Mean':(0,1),'Std':(0,4)},
                    {'name':'D','cond_type':'Normal','Mean':(0,1),'Std':(0,4)},
                    {'name':'Z','cond_type':'Normal','Mean':(0,1),'Std':(0,4)},
                    {'name':'r','cond_type':'Normal','Mean':(0,1),'Std':(0,4)}  ]                ('r')]


def run_phase_2(phase_1_path,phase_2_path,variables):

    # Problem
    e = korali.Experiment()
    e["Problem"]["Type"]  = "Hierarchical/Psi"
    e["Problem"]["Sub Problems"] = phase_1_path

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
    e["Solver"]["Type"] = "TMCMC"
    e["Solver"]["Population Size"] = 2000
    e["Solver"]["Default Burn In"] = 3;
    e["Solver"]["Target Coefficient Of Variation"] = 0.6
    e["Solver"]["Covariance Scaling"] = 0.01

    e["Console Output"]["Verbosity"] = "Detailed"
    e["File Output"]["Path"] = phase_2_path

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    k["Conduit"]["Type"] = "Concurrent"
    k["Conduit"]["Concurrent Jobs"] = 12
    k.run(e)


def run_phase_2(phase_1_path,phase_2_path):

    # ---------------------------------------------------------------------------- #
    # ---------------------------------- Setup ----------------------------------- #
    # ---------------------------------------------------------------------------- #

    e = korali.Experiment()
    e["Problem"]["Type"]  = "Hierarchical/Psi"
    e["Problem"]["Sub Problems"] = phase_1_path

    # ---------------------------------------------------------------------------- #
    # ------------------------------- Conditionals ------------------------------- #
    # ---------------------------------------------------------------------------- #
    # SIR nbin has 3 parameters: beta, gamma, [r]


    # Parameters of the conditionnal
    e["Variables"][0]["Name"] = "Psi 1" # beta mean
    e["Variables"][1]["Name"] = "Psi 2" # beta std
    e["Variables"][2]["Name"] = "Psi 3" # gamma mean
    e["Variables"][3]["Name"] = "Psi 4" # gamma std
    e["Variables"][4]["Name"] = "Psi 5" # [r] mean
    e["Variables"][5]["Name"] = "Psi 6" # [r] std

    e["Variables"][0]["Prior Distribution"] = "Uniform 0" # beta mean
    e["Variables"][1]["Prior Distribution"] = "Uniform 1" # beta std
    e["Variables"][2]["Prior Distribution"] = "Uniform 2" # gamma mean
    e["Variables"][3]["Prior Distribution"] = "Uniform 3" # gamma std
    e["Variables"][4]["Prior Distribution"] = "Uniform 4" # [r] mean
    e["Variables"][5]["Prior Distribution"] = "Uniform 5" # [r] std

    # Contidionals
    e["Problem"]["Conditional Priors"] = [ "Conditional beta", "Conditional gamma", "Conditional [r]"]

    e["Distributions"][0]["Name"] = "Conditional beta"
    e["Distributions"][0]["Type"] = "Univariate/Normal"
    e["Distributions"][0]["Mean"] = "Psi 1"
    e["Distributions"][0]["Standard Deviation"] = "Psi 2"

    e["Distributions"][1]["Name"] = "Conditional gamma"
    e["Distributions"][1]["Type"] = "Univariate/Normal"
    e["Distributions"][1]["Mean"]    = "Psi 3"
    e["Distributions"][1]["Standard Deviation"] = "Psi 4"

    e["Distributions"][2]["Name"] = "Conditional [r]"
    e["Distributions"][2]["Type"] = "Univariate/Normal"
    e["Distributions"][2]["Mean"] = "Psi 5"
    e["Distributions"][2]["Standard Deviation"] = "Psi 6"

    # ---------------------------------------------------------------------------- #
    # ---------------------------------- Priors ---------------------------------- #
    # ---------------------------------------------------------------------------- #

    e["Distributions"][3]["Name"] = "Uniform 0" # beta mean
    e["Distributions"][3]["Type"] = "Univariate/Uniform"
    e["Distributions"][3]["Minimum"] = 0.0
    e["Distributions"][3]["Maximum"] = 20.0

    e["Distributions"][4]["Name"] = "Uniform 1" # beta std
    e["Distributions"][4]["Type"] = "Univariate/Uniform"
    e["Distributions"][4]["Minimum"] = 0.0
    e["Distributions"][4]["Maximum"] = 15.0

    e["Distributions"][5]["Name"] = "Uniform 2" # gamma mean
    e["Distributions"][5]["Type"] = "Univariate/Uniform"
    e["Distributions"][5]["Minimum"] = 0.0
    e["Distributions"][5]["Maximum"] = 20.0

    e["Distributions"][6]["Name"] = "Uniform 3" # gamma std
    e["Distributions"][6]["Type"] = "Univariate/Uniform"
    e["Distributions"][6]["Minimum"] = 0.0
    e["Distributions"][6]["Maximum"] = 10.0

    e["Distributions"][7]["Name"] = "Uniform 4" # [r] mean
    e["Distributions"][7]["Type"] = "Univariate/Uniform"
    e["Distributions"][7]["Minimum"] = 0.0
    e["Distributions"][7]["Maximum"] = 5

    e["Distributions"][8]["Name"] = "Uniform 5" # [r] std
    e["Distributions"][8]["Type"] = "Univariate/Uniform"
    e["Distributions"][8]["Minimum"] = 0.0
    e["Distributions"][8]["Maximum"] = 2.5

    # ---------------------------------------------------------------------------- #
    # ---------------------------------- Solver ---------------------------------- #
    # ---------------------------------------------------------------------------- #

    e["Solver"]["Type"] = "TMCMC"
    e["Solver"]["Population Size"] = 2000
    e["Solver"]["Default Burn In"] = 3;
    e["Solver"]["Target Coefficient Of Variation"] = 0.6
    e["Solver"]["Covariance Scaling"] = 0.01

    e["Console Output"]["Verbosity"] = "Detailed"
    e["File Output"]["Path"] = phase_2_path

    # Starting Korali's Engine and running experiment
    k = korali.Engine()
    k["Conduit"]["Type"] = "Concurrent"
    k["Conduit"]["Concurrent Jobs"] = 12
    k.run(e)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='sir_int.nbin', help='Model type')
    parser.add_argument('--regions', '-r', default='cantons', help='Model type')

    args = parser.parse_args()

    model = args.model
    if args.regions == 'cantons':
        regions = CANTON_LIST
    elif args.regions == 'cantons_short':
        regions = CANTON_LIST_SHORT

    phase_1_path = ['data/'+model+'/phase_1_results/'+region+'/'+model+'/_korali_samples/' for region in regions]
    phase_2_path = 'data/'+model +'/phase_2_results/_korali_samples' 

    run_phase_2(phase_1_path,phase_2_path)



        # if var['cond_type'] == 'Normal':

        #     # Mean
        #     var_name = var['name'] + ' Mean'
        #     e["Variables"][i]["Name"] = var_name
        #     e["Variables"][i]["Prior Distribution"] = 'Uniform 'var_name

        #     j = distrib_counter+i
        #     e["Distributions"][j]["Name"] = 'Uniform '+var_name
        #     e["Distributions"][j]["Type"] = "Univariate/Uniform"
        #     e["Distributions"][j]["Minimum"] = var['Mean'][0]
        #     e["Distributions"][j]["Maximum"] = var['Mean'][1]
        #     i += 1

        #     # Std
        #     var_name = var['name'] + ' Std'

        #     e["Variables"][i]["Name"] = var_name
        #     e["Variables"][i]["Prior Distribution"] = 'Uniform 'var_name

        #     j = distrib_counter+i
        #     e["Distributions"][j]["Name"] = 'Uniform '+var_name
        #     e["Distributions"][j]["Type"] = "Univariate/Uniform"
        #     e["Distributions"][j]["Minimum"] = var['Std'][0]
        #     e["Distributions"][j]["Maximum"] = var['Std'][1]
        #     i += 1

# def add_variable(e,var_name,field,range):

#     var_name = var['name'] + ' Std'

#     e["Variables"][i]["Name"] = var_name
#     e["Variables"][i]["Prior Distribution"] = 'Uniform 'var_name

#     j = distrib_counter+i
#     e["Distributions"][j]["Name"] = 'Uniform '+var_name
#     e["Distributions"][j]["Type"] = "Univariate/Uniform"
#     e["Distributions"][j]["Minimum"] = var['Std'][0]
#     e["Distributions"][j]["Maximum"] = var['Std'][1]

#         cond_params = [ele for ele in list(var.keys()) if ele not in ['name','cond_type']] 

#         for j, cond_param in enumerate(cond_params):
#             cond_param_name = var['name']+'_'+cond_param
#             e["Variables"][i+j]["Name"] = cond_param_name
#             e["Variables"][i+1]["Prior Distribution"] = "Uniform "+cond_param_name

#     for i, cond_prior in enumerate(cond_priors):

#         e["Distributions"][0]["Name"] = cond_prior

#         e["Distributions"][0]["Type"] = "Univariate/Normal"
#         e["Distributions"][0]["Mean"] = "Psi "+str(2*i) 
#         e["Distributions"][0]["Standard Deviation"] = "Psi "+str(2*i+1)




