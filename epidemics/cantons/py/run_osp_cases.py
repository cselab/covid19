#!/usr/bin/env python3

"""Solve the ODE for many different input parameters (needed for optimal sensor placement)"""

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os
import sys
import collections

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.data.swiss_cantons import CANTON_KEYS_ALPHABETICAL, CANTON_POPULATION
from epidemics.cantons.py.model import \
        get_canton_model_data, get_canton_reference_data, \
        get_municipality_model_data, ModelData

import libepidemics

data = get_canton_model_data(include_foreign=False)
def example_run_seiin(num_days: int, inputs):
    """Runs the SEIIN model for some set of parameters and some initial conditions."""

    params = libepidemics.cantons.seiin_interventions.Parameters(
            beta =inputs[0],
            mu   =inputs[1],
            alpha=inputs[2],
            Z    =inputs[3],
            D    =inputs[4],
            theta=inputs[5],
            b0   =inputs[6],
            b1   =inputs[7],
            b2   =inputs[8],
            b3   =inputs[9])
    # Initial state.
    N0 = list(data.region_population)
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    IR0[data.key_to_index['TI']] = 1  # Ticino.
    IU0[data.key_to_index['TI']] = 10# We assume 100 unreported cases as well

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libepidemics.cantons.seiin_interventions.Solver(data.to_cpp(), verbose=False)
    return solver.solve(params, y0, num_days)




def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--days'       , type=int, default=100 , help="Number of days to evaluate.")
    parser.add_argument('--sample_file', type=str, help="Samples file; if not provided, use uniform priors")
    parser.add_argument('--samples'    , type=int, help="Samples; provide if no sample file is given")
    parser.add_argument('--no-foreign' , action='store_true',help="Disable foreign commuters from the model.")
    args = parser.parse_args(argv)

    model   = example_run_seiin
    
    cantons = 26
    days = args.days
    print ("days = ",days)

    #npar = 6 #parameters for seiir model
    npar = 10 #parameters for seiir model
    samples = 0
    if args.sample_file is None:
       print("UNIFORM PRIORS")
       samples = args.samples
       parameters = np.array([2.00,2.55,0.50,5.50,5.50,2.25, 25.0 , 60.0 , 2.00 , 2.00])
       interval   = np.array([1.50,2.45,0.50,4.50,4.50,1.75, 15.0 , 20.0 , 1.50 , 1.50])
	   
       P = np.random.uniform( 0.0, 1.0, (npar,samples))
       for s in range(samples):
           P [:,s] = (parameters-interval) + (2*interval)*P[:,s]
    
    else:
       print("NON UNIFORM PRIORS")
       sam = np.load(args.sample_file)
       #samples = sam.shape[0]
       samples = args.samples
       P = np.zeros((npar,samples))

       for s in range(samples):
           P [0:npar,s] = sam[s,0:npar]

    All_results  = np.zeros((samples,int(days),cantons))
    all_params   = np.zeros((samples,npar))
    for isim in range(samples):
        results = model(days,P[:,isim])
        for day in range(0,days):
            All_results[isim][int(day)][:] = results[day].Iu()
        all_params[isim]=P[:,isim]
        print (isim + 1,"/",samples)
    
    np.save("output_Ntheta={:05d}.npy".format(samples),All_results)
    np.save("params_Ntheta={:05d}.npy".format(samples),all_params )

if __name__ == '__main__':
    main(sys.argv[1:])
