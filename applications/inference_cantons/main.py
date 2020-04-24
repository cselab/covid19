#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import copy
from epidemics.tools.tools import import_from
import argparse

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.tools.tools import prepare_folder, save_file
from epidemics.epidemics import EpidemicsBase
import time

def repeat(v, numRegions):
  return list(np.array([v] * numRegions).T.flatten())

def mul2(v):
  v = np.array(v).astype(float)
  v[::2] *= 0.5
  return list(v)

class Models:
  SIR = "SIR"
  SEI = "SEI"

class Model( EpidemicsBase ):
  def __init__( self, data, nPropagation=100, percentages=[0.5, 0.95, 0.99], logPlot=False, **kwargs):
    self.modelName        = 'cantons'
    self.modelDescription = 'Fit SEI* with cantons on Daily Infected Data'
    self.likelihoodModel  = 'Negative Binomial'
    self.nPropagation = nPropagation
    self.percentages  = percentages
    self.logPlot = logPlot
    self.numRegions = 2

    self.model = Models.SEI
    #self.model = Models.SIR

    self.propagationData = {}

    super().__init__( **kwargs )

    self.process_data(data)

    if self.model == Models.SEI:
      self.params_fixed = {
          'R0':1.75, 'Z':2, 'D':0.8, 'tact':24}
      self.params_prior = {
          'R0':(0.5,4), 'Z':(0.1,10), 'D':(0.1,10), 'tact':(10,50)}
      self.params_to_infer = ['R0', 'Z', 'D', 'tact']
      #self.params_to_infer = []
    elif self.model == Models.SIR:
      self.params_fixed = {
          'gamma':60., 'R0':1.002}
      self.params_prior = {
          'gamma':(0.01,100), 'R0':(0.5,4)}
      self.params_to_infer = ['R0', 'gamma']
    else:
      raise NotImplementedError()


  def solve_ode(self, t_span, y0, params, t_eval):
    def smooth_trans(u0, u1, t, tact, teps):
        t0 = tact - teps
        t1 = tact + teps
        return u0 if t <= t0 else u1 if t > t1 else \
            u0 + (u1 - u0) * (1 - np.cos(np.pi/(t1 - t0)*(t - t0))) * 0.5
    def SIR(t, y):
      N = params['N']
      gamma = params['gamma']
      beta = params['R0'] * gamma
      #beta = params['beta']
      S, I = y
      c1 = beta * S * I / N
      c2 = gamma * I
      dSdt = -c1
      dIdt =  c1 - c2
      return dSdt, dIdt
    def SEI(t, y):
      N = params['N']
      Z = params['Z']
      D = params['D']
      beta = params['R0'] / D
      #beta = params['beta']
      beta = smooth_trans(beta, beta * 0.5, t, params['tact'], 0)
      S, E, I = y


      c1 = beta * S * I / N
      c2 = E / Z
      dSdt = -c1
      dEdt = c1 - c2
      dIdt = c2 - I / D
      return dSdt, dEdt, dIdt

    if self.model == Models.SEI:
      S,I = y0
      E = 0
      sol = solve_ivp(SEI, t_span=t_span, y0=[S,E,I], t_eval=t_eval)
      S,E,I = sol.y
      return [S,I]
    elif self.model == Models.SIR:
      sol = solve_ivp(SIR, t_span=t_span, y0=y0, t_eval=t_eval)
      return sol.y
    else:
      raise NotImplementedError()

  def save_data_path( self ):
      return ( self.dataFolder, self.modelName )

  def process_data(self, data):
    y = data.infected
    t = data.time
    N = data.populationSize
    I0 = y[0]
    S0 = N - I0
    y0 = S0, I0

    self.data['Model']['x-data'] = repeat(t[1:], self.numRegions)
    self.data['Model']['y-data'] = mul2(repeat(np.diff(y), self.numRegions))

    self.data['Model']['Initial Condition'] = [y0]
    self.data['Model']['Population Size'] = [N]

    T = np.ceil(t[-1])
    self.data['Propagation']['x-data'] = repeat(np.linspace(0,T,int(T+1)), self.numRegions)

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )

  def set_variables_and_distributions( self ):
    p = self.params_to_infer + ['[r]']
    for k,name in enumerate(p):
      self.e['Variables'][k]['Name'] = name
      self.e['Variables'][k]['Prior Distribution'] = 'Prior for ' + name

    self.nParameters = len(p)

    k = 0

    for name in self.params_to_infer:
      self.e['Distributions'][k]['Name'] = 'Prior for ' + name
      self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
      minmax = self.params_prior[name]
      self.e['Distributions'][k]['Minimum'] = minmax[0]
      self.e['Distributions'][k]['Maximum'] = minmax[1]
      k += 1

    self.e['Distributions'][k]['Name'] = 'Prior for [r]'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0.01
    self.e['Distributions'][k]['Maximum'] = 10.
    k += 1

  def get_params(self, p, N):
    params = {'N':N}
    for name,value in self.params_fixed.items():
      params[name] = value
    for i,name in enumerate(self.params_to_infer):
      params[name] = p[i]
    return params


  def computational_model( self, s ):
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0, = self.data['Model']['Initial Condition']
    N,  = self.data['Model']['Population Size']

    tt = [t[0]-1] + list(t)

    params = self.get_params(p, N)
    y = self.solve_ode(t_span=[0, t[-1]], y0=y0, params=params, t_eval=tt[::self.numRegions])
    S = y[0]
    Idaily = -np.diff(S)
    Idaily = mul2(repeat(Idaily, self.numRegions))

    s['Reference Evaluations'] = list(Idaily)
    s['Dispersion'] = len(Idaily) * [p[-1]]

  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0, = self.data['Model']['Initial Condition']
    N,  = self.data['Model']['Population Size']

    params = self.get_params(p, N)
    t1 = t[0::self.numRegions]
    y = self.solve_ode(t_span=[0, t[-1]], y0=y0, params=params, t_eval=t1)
    S = y[0]
    Idaily = -np.diff(S)
    Idaily = [0] + list(Idaily)
    Idaily = mul2(repeat(Idaily, self.numRegions))

    js = {}
    js['Variables'] = []

    # variables that will be visible in self.propagatedVariables
    # (does not affect sampling)
    for i in range(self.numRegions):
      js['Variables'].append({})
      js['Variables'][i]['Name']   = 'Daily Incidence {:}'.format(i)
      js['Variables'][i]['Values'] = Idaily[i::self.numRegions]

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t1)

    js['Dispersion'] = len(t1) * [p[-1]]
    s['Saved Results'] = js

  def compute_plot_intervals( self, varName, ns, ax, ylabel, cummulate=-1):
    xdata = self.data['Propagation']['x-data'][::self.numRegions]
    Ns = self.propagatedVariables[varName].shape[0]
    Nt = self.propagatedVariables[varName].shape[1]

    samples = np.zeros((Ns*ns,Nt))

    print(f"[Epidemics] Sampling from {self.likelihoodModel} for '{varName}' variable... ", end='', flush=True)

    start = time.process_time()

    if self.likelihoodModel=='Normal':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        r = self.propagatedVariables['Standard Deviation'][:,k]
        x = [ np.random.normal(m,r) for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    elif self.likelihoodModel=='Positive Normal':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        s = self.propagatedVariables['Standard Deviation'][:,k]
        t = get_truncated_normal(m,s,0,np.Inf)
        x = [ t.rvs() for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    elif self.likelihoodModel=='Negative Binomial':
      for k in range(Nt):
        m = self.propagatedVariables[varName][:,k]
        r = self.propagatedVariables['Dispersion'][:,k]
        p = p =  m/(m+r)
        x = [ np.random.negative_binomial(r,1-p) for _ in range(ns) ]
        samples[:,k] = np.asarray(x).flatten()

    else:
      sys.exit("\n[Epidemics] Likelihood not found in compute_plot_intervals.\n")

    if cummulate>0 :
      samples = np.cumsum(samples,axis=cummulate)

    elapsed = time.process_time() - start
    print(f" elapsed {elapsed:.2f} sec")

    print(f"[Epidemics] Computing quantiles... ")

    mean   = np.zeros((Nt,1))
    median = np.zeros((Nt,1))
    for k in range(Nt):
      median[k] = np.quantile( samples[:,k],0.5)
      mean[k] = np.mean( samples[:,k] )

    for p in np.sort(self.percentages)[::-1]:
      q1 = np.zeros((Nt,));
      q2 = np.zeros((Nt,));
      for k in range(Nt):
        q1[k] = np.quantile( samples[:,k],0.5-p/2)
        q2[k] = np.quantile( samples[:,k],0.5+p/2)
      ax.fill_between( xdata, q1 , q2,  alpha=0.5, label=f' {100*p:.1f}% credible interval' )

    ax.plot( xdata, mean, '-', lw=2, label='Mean', color='black')
    ax.plot( xdata, median, '--', lw=2, label='Median', color='blue')

    ax.legend(loc='upper left')
    ax.set_ylabel( ylabel )
    #ax.set_xticks( range( np.ceil( max( xdata )+1 ).astype(int) ) )
    ax.grid()
    if( self.logPlot ): ax.set_yscale('log')

    plt.draw()

    return samples

  def plot_intervals(self, region=0):
    print('[Epidemics] Compute and Plot credible intervals.')
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(self.modelDescription)

    xdata = self.data['Model']['x-data'][region::self.numRegions]
    ydata = self.data['Model']['y-data'][region::self.numRegions]

    ax  = fig.subplots( 2 )
    ax[0].plot( xdata, ydata, 'o', lw=2, label='Daily Infected(data)', color='black')

    var = 'Daily Incidence {:}'.format(region)

    self.compute_plot_intervals( var, self.nPropagation, ax[0], 'Daily Incidence' )

    #----------------------------------------------------------------------------
    z = np.cumsum(ydata)
    ax[1].plot( xdata, z, 'o', lw=2, label='Cummulative Infected (data)', color='black')

    self.compute_plot_intervals( var, self.nPropagation, ax[1], 'Cummulative number of infected', cummulate=1)

    #-----------------------------------------------------------------------------
    ax[-1].set_xlabel('time in days')

    name = "prediction{:}.png".format(str(region) if region else "")
    f = os.path.join(self.saveInfo['figures'], name);
    prepare_folder( os.path.dirname(f), clean=False )
    fig.savefig(f)

    plt.show()

    plt.close(fig)

def main():
    nSamples = 1000
    x = argparse.Namespace()
    x.dataFolder = "data/"
    x.nPropagation = 20
    x.percentages = [0.5]
    x.nThreads = 8

    data = argparse.Namespace()
    if 0: # actual data
      from epidemics.data.combined import RegionalData
      regionalData = RegionalData('switzerland')
      data.infected = regionalData.infected
      def moving_average(x, w):
        '''
        x: `numpy.ndarray`, (N)
        w: int
          Window half-width.
        Returns:
        xa: `numpy.ndarray`, (N)
          Array `x` averaged over window [-w,w].
        '''
        s = np.zeros_like(x)
        q = np.zeros_like(x)
        for i in range(len(x)):
          for j in range(max(0, i - w), min(i + w + 1, len(x))):
            if np.isfinite(x[j]):
              s[i] += x[j]
              q[i] += 1
        return s / q
      data.infected = moving_average(data.infected, 2)
      data.time = regionalData.time
      sel = len(data.infected)
      #sel = 30
      data.infected = data.infected[:sel]
      data.time = data.time[:sel]
      data.populationSize = regionalData.populationSize
      #data.populationSize = 30000
    else: # synthetic data
      N = 8e6
      I = 15
      t = range(0, 40)
      params_fixed = {
          'R0':1.75, 'Z':2, 'D':0.8, 'tact':24}
      params = {'N':N}
      params.update(params_fixed)
      m = argparse.Namespace()
      m.model = Models.SEI
      y = Model.solve_ode(m, [0,max(t)], [N-I,I], params,t)
      data.infected = y[0][0] - y[0] + I
      data.time = t
      data.populationSize = N

    x.data = data


    a = Model( **vars(x) )

    a.sample( nSamples )

    a.propagate()

    a.save()

    #a.plot_intervals()

if __name__ == "__main__":
    main()
