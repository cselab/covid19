#!/usr/bin/env python3

# Importing computational model
import sys
import os
import korali

sys.path.append('../../')
from epidemics.data.files.canton_population import CANTON_LIST, CANTON_LIST_SHORT

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
    e["Variables"][0]["Name"] = "Psi 1" # R0 mean
    e["Variables"][1]["Name"] = "Psi 2" # R0 std
    e["Variables"][2]["Name"] = "Psi 3" # gamma mean
    e["Variables"][3]["Name"] = "Psi 4" # gamma std
    e["Variables"][4]["Name"] = "Psi 5" # delta mean
    e["Variables"][5]["Name"] = "Psi 6" # delta std    
    e["Variables"][6]["Name"] = "Psi 7" # td mean
    e["Variables"][7]["Name"] = "Psi 8" # td std
    e["Variables"][6]["Name"] = "Psi 7" # [r] mean
    e["Variables"][7]["Name"] = "Psi 8" # [r] std

    e["Variables"][0]["Prior Distribution"] = "Uniform 0" # R0 mean
    e["Variables"][1]["Prior Distribution"] = "Uniform 1" # R0 std
    e["Variables"][2]["Prior Distribution"] = "Uniform 2" # gamma mean
    e["Variables"][3]["Prior Distribution"] = "Uniform 3" # gamma std
    e["Variables"][4]["Prior Distribution"] = "Uniform 4" # [r] mean
    e["Variables"][5]["Prior Distribution"] = "Uniform 5" # [r] std

    # Contidionals
    e["Problem"]["Conditional Priors"] = [ "Conditional R0", "Conditional gamma", "Conditional [r]"]

    e["Distributions"][0]["Name"] = "Conditional R0"
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

    e["Distributions"][3]["Name"] = "Uniform 0" # R0 mean
    e["Distributions"][3]["Type"] = "Univariate/Uniform"
    e["Distributions"][3]["Minimum"] = 0.0
    e["Distributions"][3]["Maximum"] = 30.0

    e["Distributions"][4]["Name"] = "Uniform 1" # R0 std
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

    model = 'sir_int.nbin'
    # cantons = ['ZH','BE','VD','GE']
    cantons = CANTON_LIST

    phase_1_path = ['data/'+model+'/phase_1_results/'+canton+'/'+model+'/_korali_samples/' for canton in cantons]
    phase_2_path = 'data/'+model +'/phase_2_results/_korali_samples' 

    run_phase_2(phase_1_path,phase_2_path)








