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
    e["Variables"][0]["Name"] = "Psi 0" # R0 mean
    e["Variables"][1]["Name"] = "Psi 1" # R0 std
    e["Variables"][2]["Name"] = "Psi 2" # gamma mean
    e["Variables"][3]["Name"] = "Psi 3" # gamma std
    e["Variables"][4]["Name"] = "Psi 4" # delta mean
    e["Variables"][5]["Name"] = "Psi 5" # delta std    
    e["Variables"][6]["Name"] = "Psi 6" # td mean
    e["Variables"][7]["Name"] = "Psi 7" # td std
    e["Variables"][8]["Name"] = "Psi 8" # [r] mean
    e["Variables"][9]["Name"] = "Psi 9" # [r] std

    e["Variables"][0]["Prior Distribution"] = "Uniform 0" # R0 mean
    e["Variables"][1]["Prior Distribution"] = "Uniform 1" # R0 std
    e["Variables"][2]["Prior Distribution"] = "Uniform 2" # gamma mean
    e["Variables"][3]["Prior Distribution"] = "Uniform 3" # gamma std
    e["Variables"][4]["Prior Distribution"] = "Uniform 4" # delta mean
    e["Variables"][5]["Prior Distribution"] = "Uniform 5" # delta std
    e["Variables"][6]["Prior Distribution"] = "Uniform 6" # td mean
    e["Variables"][7]["Prior Distribution"] = "Uniform 7" # td std    
    e["Variables"][8]["Prior Distribution"] = "Uniform 8" # [r] mean
    e["Variables"][9]["Prior Distribution"] = "Uniform 9" # [r] std
    
    # Contidionals
    e["Problem"]["Conditional Priors"] = [ "Conditional R0", "Conditional gamma", 
                                           "Conditional delta", "Conditional td",
                                           "Conditional [r]"]

    e["Distributions"][0]["Name"] = "Conditional R0"
    e["Distributions"][0]["Type"] = "Univariate/Normal"
    e["Distributions"][0]["Mean"] = "Psi 0"
    e["Distributions"][0]["Standard Deviation"] = "Psi 1"

    e["Distributions"][1]["Name"] = "Conditional gamma"
    e["Distributions"][1]["Type"] = "Univariate/Normal"
    e["Distributions"][1]["Mean"] = "Psi 2"
    e["Distributions"][1]["Standard Deviation"] = "Psi 3"

    e["Distributions"][2]["Name"] = "Conditional delta"
    e["Distributions"][2]["Type"] = "Univariate/Normal"
    e["Distributions"][2]["Mean"] = "Psi 4"
    e["Distributions"][2]["Standard Deviation"] = "Psi 5"

    e["Distributions"][3]["Name"] = "Conditional td"
    e["Distributions"][3]["Type"] = "Univariate/Normal"
    e["Distributions"][3]["Mean"] = "Psi 6"
    e["Distributions"][3]["Standard Deviation"] = "Psi 7"

    e["Distributions"][4]["Name"] = "Conditional [r]"
    e["Distributions"][4]["Type"] = "Univariate/Normal"
    e["Distributions"][4]["Mean"] = "Psi 8"
    e["Distributions"][4]["Standard Deviation"] = "Psi 9"

    # ---------------------------------------------------------------------------- #
    # ---------------------------------- Priors ---------------------------------- #
    # ---------------------------------------------------------------------------- #

    e["Distributions"][5]["Name"] = "Uniform 0" # R0 mean
    e["Distributions"][5]["Type"] = "Univariate/Uniform"
    e["Distributions"][5]["Minimum"] = 0.5
    e["Distributions"][5]["Maximum"] = 1.5

    e["Distributions"][6]["Name"] = "Uniform 1" # R0 std
    e["Distributions"][6]["Type"] = "Univariate/Uniform"
    e["Distributions"][6]["Minimum"] = 0.0
    e["Distributions"][6]["Maximum"] = 1

    e["Distributions"][7]["Name"] = "Uniform 2" # gamma mean
    e["Distributions"][7]["Type"] = "Univariate/Uniform"
    e["Distributions"][7]["Minimum"] = 0.0
    e["Distributions"][7]["Maximum"] = 0.5

    e["Distributions"][8]["Name"] = "Uniform 3" # gamma std
    e["Distributions"][8]["Type"] = "Univariate/Uniform"
    e["Distributions"][8]["Minimum"] = 0.0
    e["Distributions"][8]["Maximum"] = 0.1

    e["Distributions"][9]["Name"] = "Uniform 4" # delta mean
    e["Distributions"][9]["Type"] = "Univariate/Uniform"
    e["Distributions"][9]["Minimum"] = 0.0
    e["Distributions"][9]["Maximum"] = 1.0

    e["Distributions"][10]["Name"] = "Uniform 5" # delta std
    e["Distributions"][10]["Type"] = "Univariate/Uniform"
    e["Distributions"][10]["Minimum"] = 0.0
    e["Distributions"][10]["Maximum"] = 0.5

    e["Distributions"][11]["Name"] = "Uniform 6" # td mean
    e["Distributions"][11]["Type"] = "Univariate/Uniform"
    e["Distributions"][11]["Minimum"] = 15.
    e["Distributions"][11]["Maximum"] = 30.

    e["Distributions"][12]["Name"] = "Uniform 7" # td std
    e["Distributions"][12]["Type"] = "Univariate/Uniform"
    e["Distributions"][12]["Minimum"] = 0.0
    e["Distributions"][12]["Maximum"] = 4.0

    e["Distributions"][13]["Name"] = "Uniform 8" # [r] mean
    e["Distributions"][13]["Type"] = "Univariate/Uniform"
    e["Distributions"][13]["Minimum"] = 2.0
    e["Distributions"][13]["Maximum"] = 10.0

    e["Distributions"][14]["Name"] = "Uniform 9" # [r] std
    e["Distributions"][14]["Type"] = "Univariate/Uniform"
    e["Distributions"][14]["Minimum"] = 0.0
    e["Distributions"][14]["Maximum"] = 5.

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
    parser.add_argument('--dir', '-dir', default='./data/', help='Model type')

    args = parser.parse_args()

    model = args.model
    if args.regions == 'cantons':
        regions = CANTON_LIST
    elif args.regions == 'cantons_short':
        regions = CANTON_LIST_SHORT

    phase_1_path = [args.dir+model+'/phase_1_results/'+region+'/'+model+'/_korali_samples/' for region in regions]
    phase_2_path = args.dir+model +'/phase_2_results/_korali_samples' 

    run_phase_2(phase_1_path,phase_2_path)








