#!/usr/bin/env python3

import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import libepidemics #cpp backend

dt = 1e-3

def sird_int(y0, t_eval, N, p ):
    sird_int   = libepidemics.country.sird_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = sird_int.Solver(dp)

    params = sird_int.Parameters(R0=p[0], D=p[1], eps=p[2], tact=p[3], dtact=p[4], kbeta=p[5])
  
    initial = sird_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.1)
    
    infected   = np.zeros(len(cpp_res))
    deaths     = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S()
        deaths[idx]    = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0

    sol = np.concatenate((infected, deaths), axis=0)
    return sol


def seiird2_int(y0, t_eval, N, p ):
 
    seiird2_int = libepidemics.country.seiird2_int_reparam
    dp          = libepidemics.country.DesignParameters(N=N)
    cppsolver   = seiird2_int.Solver(dp)

    params = seiird2_int.Parameters(R0=p[0], D=p[1], Z=p[2], mu=p[3], alpha=p[4], eps=p[5], tact=p[6], dtact=p[7], kbeta=p[8])
    
    initial = seiird2_int.State(y0)
 
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
 
    infected   = np.zeros(len(cpp_res))
    deaths     = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()-entry.E()-entry.Iu()
        deaths[idx]   = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0

    sol = np.concatenate((infected, deaths), axis=0)
    return sol


def saphire_int(y0, t_eval, N, p ):
    saphire_int = libepidemics.country.saphire_int_reparam
    dp          = libepidemics.country.DesignParameters(N=N)
    cppsolver   = saphire_int.Solver(dp)

    params = saphire_int.Parameters(R0=p[0], D=p[1], Z=p[2], Y=p[3], mu=p[4], alpha=p[5], eps=p[6], tact=p[7], dtact=p[8], kbeta=p[9])
  
    initial = saphire_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
 
    infected   = np.zeros(len(cpp_res))
    deaths     = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx] = N-entry.S()-entry.E()-entry.P()-entry.Iu()
        deaths[idx]   = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0

    sol = np.concatenate((infected, deaths), axis=0)
    return sol


def seirud_int(y0, t_eval, N, p ):
    seirud_int = libepidemics.country.seirud_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = seirud_int.Solver(dp)

    params = seirud_int.Parameters(R0=p[0], D=p[1], Z=p[2], Y=p[3], alpha=p[4], eps=p[5], tact=p[6], dtact=p[7], kbeta=p[8])
  
    initial = seirud_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.01)
    
    infected        = np.zeros(len(cpp_res))
    deaths          = np.zeros(len(cpp_res))

    for idx,entry in enumerate(cpp_res):
        infected[idx]        = N-entry.S()-entry.E()-entry.P()-entry.Iu()
        deaths[idx]          = entry.D()

    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0

    sol = np.concatenate((infected, deaths), axis=0)
    return sol


def plot(total, description):
    N = len(total)

    infections = total[:int(N/2)]
    cases = np.diff(infections)
    plt.plot(infections)
    plt.ylabel('Total Infections')
    
    fname1 = "figures/{0}_cum".format(description)
    plt.savefig(fname1)
    plt.clf()
 
    plt.plot(cases)
    plt.ylabel('Daily Incidents')
    
    fname2 = "figures/{0}".format(description)
    plt.savefig(fname2)
    plt.clf()

def make_data_with_nbin_noise(total, r = 1.0):
    N = len(total)
    
    i = np.diff(total[:int(N/2)])
    d = np.diff(total[int(N/2):])
    
    pi = i/(i+r)
    xi = np.random.negative_binomial(r,1-pi)
    
    di = d/(d+r)
    xd = np.random.negative_binomial(r,1-di)

    cinfected   = np.cumsum(xi)
    cdeaths     = np.cumsum(xd)
    
    if(cinfected[0] == 0): 
        cinfected[0] = 3 # cant be 0, used to initialize infections
    sol = np.concatenate((cinfected, cdeaths), axis=0)
    
    return sol



def makefile(name, description, N, susceptible):
    fdir = "./data/{0}".format(name)
    f = open(fdir, "w+")
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
    np.random.seed(0xC0FF33)
    dispersion = 1.0

    T     = 100
    N     = 10000000
    teval = np.linspace(0, T, num=T+1)
 
    # Model Parameter
    r0    = 1.95*2
    D     = 5.2*2
    gamma = 1.0/D
    mu    = 0.9
    alpha = 0.3  # percentage reported
    Z     = 5.0  # avg latency period 
    Zl    = 3.0  # avg latency period 
    Y     = 2.0  # avg presymptomatic
    eps   = 0.02 # case fatality rate
    tact  = 20
    dtact = 7
    kbeta = 0.15/2

    beta = gamma*r0


    # Initial Conditions
    Ir0 = 10
    Iu0 = (1-alpha)/alpha*Ir0
    P0 = beta*Y*Ir0/alpha
    E0_saphire = beta*Zl*P0
    E0_seiird  = beta*Z*Ir0/alpha
    S0 = N-Ir0
    S0_saphire = N-E0_saphire-P0-Ir0-Iu0
    S0_seiird = N-E0_seiird-Ir0-Iu0
    R0 = 0
    D0 = 0

    
    # Parameter Construction
    p_sird    = [r0, D, eps, tact, dtact, kbeta]
    p_seiird  = [r0, D, Z, mu, alpha, eps, tact, dtact, kbeta]
    p_saphire = [r0, D, Zl, Y, mu, alpha, eps, tact, dtact, kbeta]
    p_seirud  = [r0, D, Zl, Y, alpha, eps, tact, dtact, kbeta]

    # SIR 
    sird_infected = sird_int((S0, Ir0, R0, D0), teval, N, p_sird)
    plot(sird_infected, "SIRD_raw")
    makefile("sird_int_raw.txt", "Synthetic SIRD Raw", N, sird_infected)

    sird_infected_rnd = make_data_with_nbin_noise(sird_infected, dispersion)
    plot(sird_infected_rnd, "SIRD_rnd")
    makefile("sird_int_rnd.txt", "Synthetic SIRD Rnd", N, sird_infected_rnd)
    
    # SEIIRD
    seiird_res = seiird2_int((S0_seiird, E0_seiird, Ir0, Iu0, R0, D0), teval, N, p_seiird)
    plot(seiird_res, "SEIIRD_raw")
    makefile("seiird2_int_raw.txt", "Synthetic SEIIRD Raw", N, seiird_res)

    seiird_res_rnd = make_data_with_nbin_noise(seiird_res, dispersion)
    plot(seiird_res_rnd, "SEIIRD_rnd")
    makefile("seiird2_int_rnd.txt", "Synthetic SEIIRD Rnd", N, seiird_res_rnd)

    # SAPHIRE 
    saphire_res = saphire_int((S0_saphire, E0_saphire, P0, Ir0, Iu0, R0, D0), teval, N, p_saphire)
    plot(saphire_res, "SAPHIRE_raw")
    makefile("saphired_int_raw.txt", "Synthetic SAPHIRE Raw", N, saphire_res)

    saphire_res_rnd = make_data_with_nbin_noise(saphire_res, dispersion)
    plot(saphire_res_rnd, "SAPHIRE_rnd")
    makefile("saphired_int_rnd.txt", "Synthetic SAPHIRE Rnd", N, saphire_res_rnd)

    # SEIRUD
    seirud_res = seirud_int((S0_saphire, E0_saphire, P0, Ir0, Iu0, R0, D0), teval, N, p_seirud)
    plot(seirud_res, "SEIRUD_raw")
    makefile("seirud_int_raw.txt", "Synthetic SEIRUD Raw", N, seirud_res)

    seirud_res_rnd = make_data_with_nbin_noise(seirud_res, dispersion)
    plot(seirud_res_rnd, "SEIRUD_rnd")
    makefile("seirud_int_rnd.txt", "Synthetic SEIRUD Rnd", N, seirud_res_rnd)


