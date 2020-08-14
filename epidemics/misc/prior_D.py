#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import scipy.stats
import scipy.integrate as integrate

# Assuming an incubation period distribution of mean 5.2 days from a separate study of early COVID-19 cases1, 
# we inferred that infectiousness started from 2.3 days (95% CI, 0.8–3.0 days) before 
# symptom onset and peaked at 0.7 days (95% CI, −0.2–2.0 days) before symptom onset (Fig. 1c). 

# The mean incubation period was 5.2 days (95% confidence interval [CI], 4.1 to 7.0), 
# with the 95th percentile of the distribution at 12.5 days.

def model_pdf(p):
  k  = p["Parameters"][0]
  th = p["Parameters"][1]

  pct40 = 14 # median of the 80 pct mild cases cases
  pct99 = 42 # duration of severe cases (3-6w)
        
  cdf40 = scipy.stats.gamma.cdf(pct40, a=k, scale=th)
  cdf99 = scipy.stats.gamma.cdf(pct99, a=k, scale=th)
 
  p["F(x)"] = - (cdf40 - .4)**2 - (cdf99 - 0.99)**2


if __name__ == '__main__':
 

    # Starting Korali's Engine
    import korali
    k = korali.Engine()

    # Creating new experiment
    e = korali.Experiment()

    # Configuring Problem
    e["Random Seed"] = 0xC0FEE
    e["Problem"]["Type"] = "Optimization"
    e["Problem"]["Objective Function"] = model_pdf

    # Defining the problem's variables.
    e["Variables"][0]["Name"] = "k"
    e["Variables"][0]["Lower Bound"] = 0.0
    e["Variables"][0]["Upper Bound"] = 50.0
 
    e["Variables"][1]["Name"] = "theta"
    e["Variables"][1]["Lower Bound"] = 0.0
    e["Variables"][1]["Upper Bound"] = 50.0

    # Configuring CMA-ES parameters
    e["Solver"]["Type"] = "Optimizer/CMAES"
    e["Solver"]["Population Size"] = 8
    e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-15
    e["Solver"]["Termination Criteria"]["Max Generations"] = 1000

    # Running Korali
    k.run(e)
