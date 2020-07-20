#!/usr/bin/env python3

import os
import numpy as np
import itertools
import pickle
import sys

sys.path.append('../../build')
import libepidemics #cpp backend

def sir_int_r0(y0, t_eval, N, p ):
    
    sir_int_r0 = libepidemics.country.sir_int_r0
    dp         = libepidemics.country.DesignParameters(N=int(N))
    cppsolver  = sir_int_r0.Solver(dp)

    params = sir_int_r0.Parameters(r0=p[0], gamma=p[1], tact=p[2], dtact=p[3], kbeta=p[4])
    
    initial = sir_int_r0.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    

    out = np.zeros((3, len(cpp_res)))

    for idx,entry in enumerate(cpp_res):
        out[0,idx] = entry.S()
        out[1,idx] = entry.I()
        out[2,idx] = entry.R()

    return out, ['S','I','R']

def seiir_int(y0, t_eval, N, p ):
    
    seiir_int = libepidemics.country.seiir_int
    dp        = libepidemics.country.DesignParameters(N=int(N))
    cppsolver = seiir_int.Solver(dp)

    params = seiir_int.Parameters(beta=p[0], mu=p[1], alpha=p[2], Z=p[3], D=p[4], tact=p[5], dtact=p[6], kbeta=p[7])
    
    initial = seiir_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    

    out = np.zeros((5, len(cpp_res)))

    for idx,entry in enumerate(cpp_res):
        out[0,idx] = entry.S()
        out[1,idx] = entry.E()
        out[1,idx] = entry.Ir()
        out[1,idx] = entry.Iu()
        out[2,idx] = entry.R()

    return out, ['S','E','Ir','Iu','R']


def produce_data(model, plist, teval, N, y0):
    data = []
    n = len(plist)
    for i,p in enumerate(plist):
        print('{}/{}'.format(i,n))
        out, fields = model(y0, teval, N, p)
        data.append((p,fields,out))

    return data


def create_dict(model, N, Y0, teval, pars, data):
    dct = {
            "Model"      : model,
            "Population" : N,
            "IC"         : Y0,
            "Teval"      : teval, 
            "Params"     : pars, 
            "Data"       : data # list of touples (params, timeseries [dim x teval] )
    }
    return dct

def store_data(filename, out):
    with open(filename, 'wb') as handle:
        pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":

    # Population Size
    N = 1e7

    # Initial Conditions
    I0 = 1
    R0 = 0
    S0 = N - I0 - R0
    
    assert(S0>0)
    Y0sir = (S0, I0, R0)             # S, I, R
    Y0seiir = (S0, 0.0, I0, 0.0, R0) # S, E, Ir, Iu, R

    # Time Window
    T = 100
    teval = np.linspace(0, T, num=T+1)

    # Adjust parameter granularity 'num'
    num = 2
    # SIR
    r0    = np.linspace( 1.0,  2.0, num=num)
    gamma = np.linspace( 0.1,  0.5, num=num)

    # SEIIR
    beta  = np.linspace( 0.1, 10.0, num=num)
    mu    = np.linspace(0.01,  1.0, num=num)
    alpha = np.linspace(0.01,  1.0, num=num)
    Z     = np.linspace( 1.0,   50, num=num)
    D     = np.linspace(1000, 8000, num=num)
 
    # Intervention
    tact  = np.linspace(10.0, 80.0, num=num)
    dtact = np.linspace( 5.0, 30.0, num=num)
    kbeta = np.linspace( 0.5,  1.0, num=num)

    r0 = [1.8]
    gamma = [0.1]
    tact = [50]
    dtact = [10]
    kbeta = [0.4]

    plist_sir   = list(itertools.product(*[r0, gamma, tact, dtact, kbeta]))
    plist_seiir = list(itertools.product(*[beta, mu, alpha, Z, D, tact, dtact, kbeta]))
    
    data_sir   = produce_data(sir_int_r0, plist_sir, teval, N, Y0sir)
    data_seiir = produce_data(seiir_int, plist_seiir, teval, N, Y0seiir)
    
    sir_dict   = create_dict("SIR with Interventions", N, Y0sir, teval, ["r0", "gamma", "tact", "dtact", "kbeta"], data_sir)
    seiir_dict = create_dict("SEIIR with Interventions", N, Y0seiir, teval, ["beta", "mu", "alpha", "Z", "D"], data_seiir)
   
    store_data("./data/SIR.pickle", sir_dict)
    store_data("./data/SEIIR.pickle", seiir_dict)
