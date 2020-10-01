#!/usr/bin/env python3

import os
import math
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import libepidemics #cpp backend
import argparse

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
    
    dt = p[6]
    w1 = math.ceil(dt)-dt
    w2 = 1.-w1

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S()
 
        if math.floor(idx+dt) < len(deaths):
            deaths[math.floor(idx+dt)] += w1 * entry.D()
        if math.ceil(idx+dt) < len(deaths):
            deaths[math.ceil(idx+dt)]  += w2 * entry.D()


    # Fix bad values
    infected[np.isnan(infected)] = 0
    deaths[np.isnan(deaths)]     = 0

    sol = np.concatenate((infected, deaths), axis=0)
    return sol


def seird_int(y0, t_eval, N, p ):
    seird_int  = libepidemics.country.seird_int_reparam
    dp         = libepidemics.country.DesignParameters(N=N)
    cppsolver  = seird_int.Solver(dp)

    params = seird_int.Parameters(R0=p[0], D=p[1], Z=p[2], eps=p[3], tact=p[4], dtact=p[5], kbeta=p[6])
  
    initial = seird_int.State(y0)
    
    cpp_res = cppsolver.solve(params, initial, t_eval=t_eval, dt = 0.1)
    
    infected   = np.zeros(len(cpp_res))
    deaths     = np.zeros(len(cpp_res))
 
    dt = p[7]
    w1 = math.ceil(dt)-dt
    w2 = 1.-w1

    for idx,entry in enumerate(cpp_res):
        infected[idx]  = N-entry.S()-entry.E()
        
        if math.floor(idx+dt) < len(deaths):
            deaths[math.floor(idx+dt)] += w1 * entry.D()
        if math.ceil(idx+dt) < len(deaths):
            deaths[math.ceil(idx+dt)]  += w2 * entry.D()


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

    dt = p[9]
    w1 = math.ceil(dt)-dt
    w2 = 1.-w1

    for idx,entry in enumerate(cpp_res):
        infected[idx] = (N-entry.S()-entry.E())*p[4]
 
        if math.floor(idx+dt) < len(deaths):
            deaths[math.floor(idx+dt)] += w1 * entry.D()
        if math.ceil(idx+dt) < len(deaths):
            deaths[math.ceil(idx+dt)]  += w2 * entry.D()


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

    dt = p[10]
    w1 = math.ceil(dt)-dt
    w2 = 1.-w1

    for idx,entry in enumerate(cpp_res):
        infected[idx] = (N-entry.S()-entry.E()-entry.P())*p[5]
 
        if math.floor(idx+dt) < len(deaths):
            deaths[math.floor(idx+dt)] += w1 * entry.D()
        if math.ceil(idx+dt) < len(deaths):
            deaths[math.ceil(idx+dt)]  += w2 * entry.D()


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
    infectedu       = np.zeros(len(cpp_res))
    deaths          = np.zeros(len(cpp_res))

    dt = p[9]
    w1 = math.ceil(dt)-dt
    w2 = 1.-w1

    for idx,entry in enumerate(cpp_res):
        infected[idx]        = (N-entry.S()-entry.E()-entry.P())*p[4]
 
        if math.floor(idx+dt) < len(deaths):
            deaths[math.floor(idx+dt)] += w1 * entry.D()
        if math.ceil(idx+dt) < len(deaths):
            deaths[math.ceil(idx+dt)]  += w2 * entry.D()


    # Fix bad values
    infected[np.isnan(infected)]   = 0
    deaths[np.isnan(deaths)]       = 0  
    
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
    ci = total[:int(N/2)]
    cd = total[int(N/2):]
    
    i = np.diff(ci)
    d = np.diff(cd)
    
    pi = i/(i+r)
    xi = np.random.negative_binomial(r,1-pi)
    
    di = d/(d+r)
    xd = np.random.negative_binomial(r,1-di)

    cinfected    = ci[0] + np.cumsum(xi)
    cinfected[0] = ci[0]

    cdeaths    = cd[0] + np.cumsum(xd)
    cdeaths[0] = cd[0]
    
    if(cinfected[0] == 0): 
        print("WARNING, CINFECTED[0] == 0")
        cinfected[0] = 1 # cant be 0, used to initialize infections
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

    parser = argparse.ArgumentParser(description='COVID project.')
    parser.add_argument('-I0', type=int, default=10, help='Infections IC')
    parser.add_argument('-D', type=float, default=5.2, help='Infectious period duration')
    parser.add_argument('-mu', type=float, default=0.6, help='Transmission reduction factor')
    parser.add_argument('-alpha', type=float, default=0.3, help='Reporting rate')
    parser.add_argument('-Z', type=float, default=5.2, help='Latency Period before Presymptomatic')
    parser.add_argument('-Zl', type=float, default=2.9, help='Latency Period before Presymptomatic')
    parser.add_argument('-Y', type=float, default=2.28, help='Presymptomatic Period')
    parser.add_argument('-eps', type=float, default=0.05, help='Case fatality rate')
    parser.add_argument('-delay', type=float, default=7, help='Delayed death')
    parser.add_argument('-tact', type=float, default=20, help='Intervention time')
    parser.add_argument('-dtact', type=float, default=14, help='Intervention duration')
    parser.add_argument('-kbeta', type=float, default=0.2, help='Intervention reduction factor')
    parser.add_argument('-dispersion', type=float, default=50, help='NBIN dispersion')
    parser.add_argument('-r0sir', type=float, default=2.0, help='R0')
    parser.add_argument('-r0seir', type=float, default=3.0, help='R0')
    parser.add_argument('-r0seiir', type=float, default=5.0, help='R0')
    parser.add_argument('-r0seiru', type=float, default=1.75, help='R0')
    parser.add_argument('-r0saphire', type=float, default=1.15, help='R0')
    parser.add_argument('-seed', type=int, default=1337, help='Random Seed')
    args = parser.parse_args()
    
    np.random.seed(args.seed)
    dispersion = args.dispersion

    T     = 100
    N     = 10000000
    teval = np.linspace(0, T, num=T+1)
 
    # Model Parameter
    r0_sir     = args.r0sir
    r0_seir    = args.r0seir
    r0_seiir   = args.r0seiir
    r0_seiru   = args.r0seiru
    r0_saphire = args.r0saphire

    D     = args.D
    mu    = args.mu
    alpha = args.alpha  
    Z     = args.Z  
    Zl    = args.Zl
    Y     = args.Y
    eps   = args.eps
    delay = args.delay
    tact  = args.tact
    dtact = args.dtact
    kbeta = args.kbeta

    # Initial Conditions
    Ir0 = args.I0
    Iu0 = (1-alpha)/alpha*Ir0
    P0_seiru   = r0_seiru/D*Y*Ir0/alpha
    P0_saphire = r0_saphire/D*Y*Ir0/alpha
    E0_seir    = r0_seir/D*Z*Ir0
    E0_seiir   = r0_seiir/D*Z*Ir0/alpha
    E0_saphire = r0_saphire/D*Zl*P0_saphire
    E0_seiru   = r0_seiru/D*Zl*P0_seiru
    S0 = N-Ir0
    S0_seir  = N-E0_seir-Ir0
    S0_seiru = N-E0_seiru-P0_seiru-Ir0-Iu0
    S0_saphire = N-E0_saphire-P0_saphire-Ir0-Iu0
    S0_seiir = N-E0_seiir-Ir0-Iu0
    R0 = 0
    D0 = 0
    
    # Parameter Construction
    p_sird    = [r0_sir, D, eps, tact, dtact, kbeta, delay]
    p_seird   = [r0_seir, D, Z, eps, tact, dtact, kbeta, delay]
    p_seiird  = [r0_seiir, D, Z, mu, alpha, eps, tact, dtact, kbeta, delay]
    p_seirud  = [r0_seiru, D, Zl, Y, alpha, eps, tact, dtact, kbeta, delay]
    p_saphire = [r0_saphire, D, Zl, Y, mu, alpha, eps, tact, dtact, kbeta, delay]

    # SIRD
    sird_infected = sird_int((S0, Ir0, R0, D0), teval, N, p_sird)
    plot(sird_infected, "SIRD_raw")
    makefile("sirdelay_int_raw.txt", "Synthetic SIRD Raw", N, sird_infected)

    sird_infected_rnd = make_data_with_nbin_noise(sird_infected, dispersion)
    plot(sird_infected_rnd, "SIRD_rnd")
    makefile("sirdelay_int_rnd.txt", "Synthetic SIRD Rnd", N, sird_infected_rnd)
 
    # SEIRD
    seird_infected = seird_int((S0_seir, E0_seir, Ir0, R0, D0), teval, N, p_seird)
    plot(seird_infected, "SEIRD_raw")
    makefile("seirdelay_int_raw.txt", "Synthetic SEIRD Raw", N, seird_infected)

    seird_infected_rnd = make_data_with_nbin_noise(seird_infected, dispersion)
    plot(seird_infected_rnd, "SEIRD_rnd")
    makefile("seirdelay_int_rnd.txt", "Synthetic SEIRD Rnd", N, seird_infected_rnd)
    
    # SEIIRD
    seiird_res = seiird2_int((S0_seiir, E0_seiir, Ir0, Iu0, R0, D0, 0., 0.), teval, N, p_seiird)
    plot(seiird_res, "SEIIRD_raw")
    makefile("seiirdelay_int_raw.txt", "Synthetic SEIIRD Raw", N, seiird_res)

    seiird_res_rnd = make_data_with_nbin_noise(seiird_res, dispersion)
    plot(seiird_res_rnd, "SEIIRD_rnd")
    makefile("seiirdelay_int_rnd.txt", "Synthetic SEIIRD Rnd", N, seiird_res_rnd)

    # SAPHIRE 
    saphire_res = saphire_int((S0_saphire, E0_saphire, P0_saphire, Ir0, Iu0, R0, D0, 0., 0.), teval, N, p_saphire)
    plot(saphire_res, "SAPHIRE_raw")
    makefile("saphiredelay_int_raw.txt", "Synthetic SAPHIRE Raw", N, saphire_res)

    saphire_res_rnd = make_data_with_nbin_noise(saphire_res, dispersion)
    plot(saphire_res_rnd, "SAPHIRE_rnd")
    makefile("saphiredelay_int_rnd.txt", "Synthetic SAPHIRE Rnd", N, saphire_res_rnd)

    # SEIRUD
    seirud_res = seirud_int((S0_seiru, E0_seiru, P0_seiru, Ir0, Iu0, R0, D0), teval, N, p_seirud)
    plot(seirud_res, "SEIRUD_raw")
    makefile("seirudelay_int_raw.txt", "Synthetic SEIRUD Raw", N, seirud_res)

    seirud_res_rnd = make_data_with_nbin_noise(seirud_res, dispersion)
    plot(seirud_res_rnd, "SEIRUD_rnd")
    makefile("seirudelay_int_rnd.txt", "Synthetic SEIRUD Rnd", N, seirud_res_rnd)
