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
  scale  = p["Parameters"][0]

  mu    = 5.2
  pct95 = 12.5 # duration of severe cases (3-6w)
    
  mean, var = scipy.stats.weibull_min.stats(c=shape, scale = scale)
  cdf05 = 1-np.exp(-scale*pct05**shape)
  cdf95 = 1-np.exp(-scale*pct95**shape)
 
  print(mean)
  print(cdf95,flush=True)
  p["F(x)"] = - (mean - mu)**2 - (cdf95 - 0.95)**2


if __name__ == '__main__':
 

    # Starting Korali's Engine
    import korali
    k = korali.Engine()

    # Creating new experiment
    e = korali.Experiment()

    # Configuring Problem
    #e["Random Seed"] = 0xC0FEE
    e["Problem"]["Type"] = "Optimization"
    e["Problem"]["Objective Function"] = model_pdf

    # Defining the problem's variables.
 
    e["Variables"][0]["Name"] = "scale"
    e["Variables"][0]["Lower Bound"] = 0.0
    e["Variables"][0]["Upper Bound"] = 5.0

    # Configuring CMA-ES parameters
    e["Solver"]["Type"] = "Optimizer/CMAES"
    e["Solver"]["Population Size"] = 16
    e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-15
    e["Solver"]["Termination Criteria"]["Max Generations"] = 1000

    # Running Korali
    k.run(e)
