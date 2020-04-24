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

import libsolver

class Level:
    canton = "canton"
    municipality = "municipality"

def example_run_seiin(data: ModelData, num_days: int, inputs):
    """Runs the SEIIN model for some set of parameters and some initial conditions."""

    # Parameters.
    # Li2020, Table 1
    params = libsolver.solvers.seiin.Parameters(
            beta=inputs[0], mu=inputs[1], alpha=inputs[2], Z=inputs[3], D=inputs[4], theta=inputs[5])
            #beta=1.12, mu=0., alpha=1., Z=3.69, D=3.47, theta=1.36)

    # Initial state.
    N0 = list(data.region_population)
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    if 'TI' in data.key_to_index:
        IR0[data.key_to_index['TI']] = 1  # Ticino.
        IU0[data.key_to_index['TI']] = 100# We assume 100 unreported cases as well
    else:
        IR0[data.key_to_index['MUN-5192']] = 1  # Lugano.
        IU0[data.key_to_index['MUN-5192']] = 100# We assume 100 unreported cases as well

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libsolver.solvers.seiin.Solver(data.to_cpp(), verbose=False)
    return solver.solve(params, y0, num_days)


def example_run_seii_c(data, num_days,inputs):
    """Runs the SEII_C model for some set of parameters and some initial conditions."""
    # Parameters. Note that nu is relative to beta.
    params = libsolver.solvers.seii_c.Parameters(
            beta=inputs[0], mu=inputs[1], alpha=inputs[2], Z=inputs[3], D=inputs[4])
            #beta=3.0, nu=1.0, alpha=0.6, Z=3.69, D=3.47)

    # Initial state.
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    # if 'TI' in data.key_to_index:
    #     IU0[data.key_to_index['TI']] = 10  # Ticino.
    # else:
    #     IU0[data.key_to_index['MUN-5192']] = 10  # Lugano.

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(data.region_population, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0

    # Run the ODE solver.
    solver = libsolver.solvers.seii_c.Solver(data.to_cpp(), verbose=False)
    return solver.solve(params, y0, num_days)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('days', type=int, default=50, help="Number of days to evaluate.")
    parser.add_argument('samples', type=int, default=100, help="Number of Monte Carlo samples.")
    parser.add_argument('--TMCMC', type=int,  help="Use TMCMC samples or not.")
    parser.add_argument('--no-foreign', action='store_true', help="Disable foreign commuters from the model.")
    parser.add_argument('--level', type=str, choices=(Level.canton, Level.municipality), default='canton', help="Level of details.")
    parser.add_argument('--model', type=str, choices=('seiin', 'seii_c'), default='seiin', help="Model.")
    args = parser.parse_args(argv)

    if args.level == Level.canton:
        model_data = get_canton_model_data(include_foreign=not args.no_foreign)
        ref_data = get_canton_reference_data()
    elif args.level == Level.municipality:
        model_data = get_municipality_model_data()
        ref_data = None
    else:
        assert False

    if args.model == 'seiin':
        model = example_run_seiin
        samples = args.samples
        parameters = np.array([1.50,1.00,0.50,3.69,3.47,1.36, 0.3, 0.1])
        interval   = np.array([0.50,0.50,0.10,0.30,0.30,0.20, 0.0, 0.0])
    else:
        model_data.Mij *= 0.0
        model = example_run_seii_c
        parameters = np.array([3.00,1.00,0.60,3.69,3.47, 0.3  ,0.1])
        interval   = np.array([0.50,1.00,0.20,0.00,0.50, 0.0  ,0.0])


    assert args.level == Level.canton

    cantons = 26
    c = 1
    days = args.days
    print ("days = ",days)
    print ("samples = ",args.samples)

    npar = len(parameters)
    
    if args.TMCMC == 0:
       print("UNIFORM PRIORS")
       P = np.random.uniform( 0.0, 1.0, (npar,samples))
       for s in range(samples):
           P [:,s] = (parameters-interval) + (2*interval)*P[:,s]
    else:
       print("NON UNIFORM PRIORS")
       P = np.zeros((npar,samples))
       sam = np.load("samples.npy")
       for s in range(samples):
           P [0:6,s] = sam[s,:]
       


    All_results  = np.zeros((samples,int(days/c),cantons))
    reported     = np.zeros((samples,int(days/c),cantons))
    all_params   = np.zeros((samples,npar))
    for isim in range(samples):
        results = model(model_data,days,P[:,isim])
        for day in range(0,days,c):
            All_results[isim][int(day/c)][:] = results[day].Iu()
            reported   [isim][int(day/c)][:] = results[day].E()
        reported[isim,:,:] *= (P[2,isim] /P[3,isim])
        all_params[isim]=P[:,isim]
        print (isim + 1,"/",samples)
    
    np.save("output_Ntheta={:05d}.npy".format(samples)   ,All_results)
    np.save("params_Ntheta={:05d}.npy".format(samples)   ,all_params )
    np.save("reported.npy" ,reported   )

if __name__ == '__main__':
    main(sys.argv[1:])
