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
from epidemics.data.combined import RegionalData


class Model( EpidemicsBase ):
  def __init__( self, **kwargs ):
    self.modelName        = 'cantons'
    self.modelDescription = 'Fit SEI* with cantons on Daily Infected Data'
    self.likelihoodModel  = 'Negative Binomial'

    self.futureDays   = kwargs.pop('futureDays', 2)
    self.nPropagation = kwargs.pop('nPropagation', 100)
    self.logPlot      = kwargs.pop('logPlot', False)
    self.percentages  = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess   = kwargs.pop('preprocess', False)

    self.regionalData = RegionalData( 'switzerland',self.preprocess)
    self.propagationData = {}

    super().__init__( **kwargs )

    self.process_data()

    self.params_fixed = {"beta":80., "gamma":80.}
    self.params_prior = {"beta":(1,100), "gamma":(1,100)}
    #self.params_to_infer = ["beta", "gamma"]
    self.params_to_infer = ["beta"]

  def solve_model(self, t_span, y0, params, t_eval):
    def rhs(t, y):
      N = params['N']
      beta = params['beta']
      gamma = params['gamma']
      S, I = y
      c1 = beta * S * I / N
      c2 = gamma * I
      dSdt = -c1
      dIdt =  c1 - c2
      return dSdt, dIdt
    sol = solve_ivp(rhs, t_span=t_span, y0=y0, t_eval=t_eval)
    return sol.y

  def save_data_path( self ):
      return ( self.dataFolder, self.modelName )

  def process_data( self ):
    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize
    I0 = y[0]
    S0 = N - I0
    y0 = S0, I0

    self.data['Model']['x-data'] = t[1:]
    self.data['Model']['y-data'] = np.diff( y[0:])

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))

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
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()

    params = self.get_params(p, N)
    y = self.solve_model(t_span=[0, t[-1]], y0=y0, params=params, t_eval=tt)
    S = y[0]
    Idaily = -np.diff(S)

    s['Reference Evaluations'] = list(Idaily)
    s['Dispersion'] = len(Idaily) * [p[-1]]

  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    params = self.get_params(p, N)
    y = self.solve_model(t_span=[0, t[-1]], y0=y0, params=params, t_eval=t)
    S = y[0]
    Idaily = -np.diff(S)
    Idaily = [0] + list(Idaily)

    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][0]['Name']   = 'Daily Incidence'
    js['Variables'][0]['Values'] = Idaily

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    js['Dispersion'] = len(Idaily) * [p[-1]]

    s['Saved Results'] = js

  def plot_intervals( self ):
    print('[Epidemics] Compute and Plot credible intervals.')
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(self.modelDescription)

    ax  = fig.subplots( 2 )
    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', 20, ax[0], 'Daily Incidence' )

    #----------------------------------------------------------------------------
    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected (data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', 20, ax[1], 'Cummulative number of infected', cummulate=1)

    #-----------------------------------------------------------------------------
    ax[-1].set_xlabel('time in days')

    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    plt.show()

    plt.close(fig)

def main():
    nSamples = 500
    x = argparse.Namespace()
    x.dataFolder = "data/"
    x.nPropagation = 100
    x.percentages = [0.5]
    x.nThreads = 4

    a = Model( **vars(x) )

    a.sample( nSamples )

    a.propagate()

    a.save()

    #a.plot_intervals()

if __name__ == "__main__":
    main()
