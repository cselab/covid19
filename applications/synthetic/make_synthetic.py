#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

import libepidemics #cpp backend

def sir_beta(y0, t_eval, N, p ):
    
    sir       = libepidemics.country.sir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = sir.Solver(data)

    params = sir.Parameters(beta=p[0], gamma=p[1])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    initial = sir.State(y0cpp)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    Svec = np.zeros(len(cpp_res)+1)
    Svec[0] = N
    for idx,entry in enumerate(cpp_res):
        Svec[idx+1] = entry.S()

    # get newly infected
    return Svec


def sir_int_r0(y0, t_eval, N, p ):
    
    sir_int_r0 = libepidemics.country.sir_int_r0
    data       = libepidemics.country.ModelData(N=N)
    cppsolver  = sir_int_r0.Solver(data)

    params = sir_int_r0.Parameters(r0=p[0], gamma=p[1], tact=p[2], dtact=p[3], kbeta=p[4])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    initial = sir_int_r0.State(y0cpp)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    Svec = np.zeros(len(cpp_res)+1)

    Svec[0] = N
    for idx,entry in enumerate(cpp_res):
        Svec[idx+1] = entry.S()

    # get newly infected
    return Svec


def make_data_with_mul_nrm_noise(infected, sig = 1.0):
    cases = np.diff(infected)
    # Add multiplicative noise, round to integers and remove negative values
    cstar = cases + cases*np.random.normal(0, sig, len(cases))
    cstar = np.round(cstar)
    cstar[cstar < 1] = 1
    Sstar = np.cumsum(cstar)
    return Sstar

 
def plot(total):
    cases = np.diff(total)
    plt.plot(total)
    plt.ylabel('Total Infections')
    plt.show()
 
    plt.plot(cases)
    plt.ylabel('Daily Incidents')
    plt.show()


def makefile(name, description, N, susceptible):
    f = open(name, "a")
    f.write(description)
    f.write(os.linesep)
    f.write(str(N))
    f.write(os.linesep)
    f.write(str(len(susceptible)))
    f.write(os.linesep)
    for s in susceptible:
        f.write(str(s))
        f.write(os.linesep)
    f.close()

if __name__ == "__main__":

    np.random.seed(1337)

    N = 10000000

    I0 = 1
    S0 = N - I0
    T  = 90

    teval = np.linspace(0, T, num=T+1)

    r0    = 1.6
    gamma = 0.3
    tact  = 40
    dtact = 10
    kbeta = 0.5

    beta = gamma*r0
    psir = [beta, gamma]

    sirS     = sir_beta((S0, I0), teval, N, psir)
    print(sirS)
    sirCases = -np.diff(sirS)
    print(sirCases)
    sirInfected = np.cumsum(sirCases)
    print(sirInfected)
    makefile("sir_raw.txt", "Synthetic SIR Raw", N, sirInfected)
    
    p = [r0, gamma, tact, dtact, kbeta]
    
    S        = sir_int_r0((S0, I0), teval, N, p)
    cases    = -np.diff(S)
    infected = np.cumsum(cases)
    
    print("Plotting unprocessed solver output..")
    plot(infected)
    
    infectedRand = make_data_with_mul_nrm_noise(infected, 1)
    print("Plotting randomized solver output..")
    plot(infectedRand)

    makefile("sir_int_r0_raw.txt", "Synthetic Raw", N, infected)
    makefile("sir_int_r0_rndm.txt", "Synthetic Rnd", N, infectedRand)
