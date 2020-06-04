#!/usr/bin/env python3

import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import libepidemics #cpp backend


def sir_beta(y0, t_eval, N, p ):
    
    sir       = libepidemics.country.sir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = sir.Solver(data)

    params = sir.Parameters(beta=p[0], gamma=p[1])
    
    initial = sir.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected = np.zeros(len(cpp_res))
    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()
    
    return infected


def sir_int_r0(y0, t_eval, N, p ):
    
    sir_int_r0 = libepidemics.country.sir_int_r0
    data       = libepidemics.country.ModelData(N=N)
    cppsolver  = sir_int_r0.Solver(data)

    params = sir_int_r0.Parameters(r0=p[0], gamma=p[1], tact=p[2], dtact=p[3], kbeta=p[4])
    
    initial = sir_int_r0.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected = np.zeros(len(cpp_res))
    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()
 
    return infected


def seir(y0, t_eval, N, p ):
    seir      = libepidemics.country.seir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seir.Solver(data)

    params = seir.Parameters(beta=p[0], gamma=p[1], a=p[2])
    
    initial = seir.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected = np.zeros(len(cpp_res))
    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()-entry.E()
    
    return infected


def seir_int(y0, t_eval, N, p ):
    seir_int  = libepidemics.country.seir_int
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seir_int.Solver(data)

    params = seir_int.Parameters(beta=p[0], gamma=p[1], a=p[2], tact=p[3], dtact=p[4], kbeta=p[5])
    
    initial = seir_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected = np.zeros(len(cpp_res))
    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()-entry.E()
    
    return infected


def seiir(y0, t_eval, N, p ):
    seiir     = libepidemics.country.seiir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seiir.Solver(data)

    params = seiir.Parameters(beta=p[0], mu=p[1], alpha=p[2], Z=p[3], D=p[4])
    
    initial = seiir.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected = np.zeros(len(cpp_res))
    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()-entry.E()-entry.Iu()
    
    return infected


def seiir_int(y0, t_eval, N, p ):
    seiir_int = libepidemics.country.seiir_int
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seiir_int.Solver(data)

    params = seiir_int.Parameters(beta=p[0], mu=p[1], alpha=p[2], Z=p[3], D=p[4], tact=p[5], dtact=p[6], kbeta=p[7])
    
    initial = seiir_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected = np.zeros(len(cpp_res))
    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()-entry.E()-entry.Iu()
    
    return infected


def make_data_with_mul_nrm_noise(infected, sig = 1.0):
    cases = np.diff(infected)
    # Add multiplicative noise, round to integers and remove negative values
    cstar = cases + cases*np.random.normal(0, sig, len(cases))
    cstar = np.round(cstar)
    cstar[cstar < 1] = 1
    Sstar = np.cumsum(cstar)
    return Sstar

 
def plot(total, description):
    cases = np.diff(total)
    plt.plot(total)
    plt.ylabel('Total Infections')
    
    fname1 = "figures/{0}_cum".format(description)
    plt.savefig(fname1)
    plt.clf()
 
    plt.plot(cases)
    plt.ylabel('Daily Incidents')
    
    fname2 = "figures/{0}".format(description)
    plt.savefig(fname2)
    plt.clf()


def makefile(name, description, N, susceptible):
    f = open(name, "w+")
    f.seek(0) # replace old content
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

    # 'global' params
    np.random.seed(1337)
    noise = 0.1

    T     = 90
    N     = 10000000
    teval = np.linspace(0, T, num=T+1)

    # Initial Conditions
    Ir0 = 1
    Iu0 = 0
    S0 = N - Ir0
    E0 = 0
    R0 = 0
    
    # Model Parameter
    r0    = 1.75
    gamma = 1.0/5.2
    mu    = 1.0/2.0 # inverse is avg incubation period
    alpha = 0.7 # percentage reported
    Z     = 0.1 # avg latency period 
    D     = 10.0 # avg duration of infection
    tact  = 40
    dtact = 10
    kbeta = 0.3

    beta = gamma*r0 # 0.3365
    
    # Parameter Construction
    p_sir       = [beta, gamma]
    p_sir_r0    = [r0, gamma]
    p_sir_int   = [r0, gamma, tact, dtact, kbeta]
    p_seir      = [beta, gamma, mu]
    p_seir_int  = [beta, gamma, mu, tact, dtact, kbeta]
    p_seiir     = [beta, mu, alpha, Z, D]
    p_seiir_int = [beta, mu, alpha, Z, D, tact, dtact, kbeta]

    # SIR 
    sir_infected = sir_beta((S0, Ir0, R0), teval, N, p_sir)
    plot(sir_infected, "SIR")
    makefile("sir_raw.txt", "Synthetic SIR Raw", N, sir_infected)

    sir_infected_rnd = make_data_with_mul_nrm_noise(sir_infected, noise)
    plot(sir_infected_rnd, "SIR_rnd")
    makefile("sir_rnd.txt", "Synthetic SIR Rnd", N, sir_infected_rnd)
 
    # SIR with interventions
    sir_infected_int = sir_int_r0((S0, Ir0, R0), teval, N, p_sir_int)
    plot(sir_infected_int, "SIR_int")
    makefile("sir_int_r0_raw.txt", "Synthetic SIR with Interventions Raw", N, sir_infected_int)
    
    sir_infected_int_rnd = make_data_with_mul_nrm_noise(sir_infected_int, noise)
    plot(sir_infected_int_rnd,"SIR_rnd")
    makefile("sir_int_r0_rnd.txt", "Synthetic SIR with Interventions Rnd", N, sir_infected_int_rnd)

    # SEIR
    seir_infected = seir((S0, E0, Ir0, R0), teval, N, p_seir)
    plot(seir_infected, "SEIR_raw")
    makefile("seir_raw.txt", "Synthetic SEIR Raw", N, seir_infected)
 
    seir_infected_rnd = make_data_with_mul_nrm_noise(sir_infected, noise)
    plot(seir_infected_rnd,"SEIR_rnd")
    makefile("seir_rnd.txt", "Synthetic SEIR with Interventions Rnd", N, seir_infected_rnd)

    # SEIR with Interventions
    seir_infected_int = seir_int((S0, E0, Ir0, R0), teval, N, p_seir_int)
    plot(seir_infected_int, "SEIR_raw")
    makefile("seir_int_raw.txt", "Synthetic SEIR with Interventions Raw", N, seir_infected_int)

    seir_infected_int_rnd = make_data_with_mul_nrm_noise(seir_infected_int, noise)
    plot(seir_infected_int_rnd,"SEIR_int_rnd")
    makefile("seir_int_rnd.txt", "Synthetic SEIR with Interventions Rnd", N, seir_infected_int_rnd)

    # SEIIR 
    seiir_infected = seiir((S0, E0, Ir0, Iu0, R0), teval, N, p_seiir)
    plot(seiir_infected, "SEIIR_raw")
    makefile("seiir_raw.txt", "Synthetic SEIIR Raw", N, seiir_infected)
 
    seiir_infected_rnd = make_data_with_mul_nrm_noise(seiir_infected, noise)
    plot(seiir_infected_rnd,"SEIIR_rnd")
    makefile("seiir_rnd.txt", "Synthetic SEIR Rnd", N, seiir_infected_rnd)


    # SEIIR with Interventions
    seiir_infected_int = seiir_int((S0, E0, Ir0, Iu0, R0), teval, N, p_seiir_int)
    plot(seiir_infected_int, "SEIIR_int_raw")
    makefile("seiir_int_raw.txt", "Synthetic SEIIR with Interventions Raw", N, seiir_infected_int)

    seiir_infected_int_rnd = make_data_with_mul_nrm_noise(seiir_infected_int, noise)
    plot(seiir_infected_int_rnd,"SEIIR_int_rnd")
    makefile("seiir_int_rnd.txt", "Synthetic SEIIR with Interventions Rnd", N, seiir_infected_int_rnd)
