#!/usr/bin/env python3

import os
import numpy as np
import itertools
import pickle

import libepidemics #cpp backend


def sir_int_r0(y0, t_eval, N, p ):
    
    sir_int_r0 = libepidemics.country.sir_int_r0
    data       = libepidemics.country.ModelData(N=int(N))
    cppsolver  = sir_int_r0.Solver(data)

    params = sir_int_r0.Parameters(r0=p[0], gamma=p[1], tact=p[2], dtact=p[3], kbeta=p[4])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    initial = sir_int_r0.State(y0cpp)
    
    cpp_res = cppsolver.solve_params_ad(params, initial, t_eval=t_eval, dt = 0.01)
    

    out = np.zeros((3, len(cpp_res)))

    for idx,entry in enumerate(cpp_res):
        out[0,idx] = entry.S().val()
        out[1,idx] = entry.I().val()
        out[2,idx] = entry.R().val()

    return out

def produce_params(r0, gamma, tact, dtact, kbeta):
    ps = [r0, gamma, tact, dtact, kbeta]
    paramslst = list(itertools.product(*ps))
    return paramslst

def produce_data(model, plist, teval, N, y0):
    data = []
    for p in plist:
        out = model(y0, teval, N, p)
        data.append((p,out))

    return data

def store_data(filename, out):
    with open(filename, 'wb') as handle:
        pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":

    # Population Size
    N = 1e7

    # Initial Conditions
    I0 = 1
    S0 = N - I0 
    
    assert(S0>0)
    Y0 = (S0, I0)

    # Time Window
    T = 100
    teval = np.linspace(0, T, num=T+1)

    r0    = np.linspace( 1.0,  2.0, num=2)
    gamma = np.linspace( 0.1,  0.5, num=2)
    tact  = np.linspace(10.0, 80.0, num=2)
    dtact = np.linspace( 5.0, 30.0, num=2)
    kbeta = np.linspace( 0.5,  1.0, num=2)

    plist = produce_params(r0, gamma, tact, dtact, kbeta)
    data  = produce_data(sir_int_r0, plist, teval, N, Y0)

    output = {
            "Model"      : "SIR with Interventions",
            "Population" : N,
            "IC"         : Y0,
            "Teval"      : teval,
            "Params"     : ["r0", "gamma" , "tact", "dtact", "kbeta"],
            "Data"       : data
    }
    
    store_data("SIR.pickle", output)
