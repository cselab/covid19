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

    self.country      = kwargs.pop('country', 'switzerland')
    self.futureDays   = kwargs.pop('futureDays', 2)
    self.nPropagation = kwargs.pop('nPropagation', 100)
    self.logPlot      = kwargs.pop('logPlot', False)
    self.nValidation  = kwargs.pop('nValidation', 0)
    self.percentages  = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess   = kwargs.pop('preprocess', False)

    self.regionalData = RegionalData( self.country,self.preprocess)
    self.propagationData = {}

    super().__init__( **kwargs )

    self.process_data()

  def solve_model(self, t_span, y0, N, p, t_eval):
    def sir_rhs(t, y, N, p):
      S, I = y
      beta, gamma = p[:2]
      c1 = beta * S * I / N
      c2 = gamma * I
      dSdt = -c1
      dIdt =  c1 - c2
      return dSdt, dIdt
    sol = solve_ivp(sir_rhs, t_span=t_span, y0=y0, args=(N,p), t_eval=t_eval)
    return sol.y

  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )

  def process_data( self ):
    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize
    I0 = y[0]
    S0 = N - I0
    y0 = S0, I0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = np.diff( y[0:])
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = np.diff( y[0:-self.nValidation] )
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = np.diff( y[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )

  def set_variables_and_distributions( self ):
    p = ['beta','gamma','[r]']
    for k,x in enumerate(p):
      self.e['Variables'][k]['Name'] = x
      self.e['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

    self.nParameters = len(p)

    k=0
    self.e['Distributions'][k]['Name'] = 'Prior for beta'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 1.
    self.e['Distributions'][k]['Maximum'] = 100.
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for gamma'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 1.
    self.e['Distributions'][k]['Maximum'] = 100.
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for [r]'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0.01
    self.e['Distributions'][k]['Maximum'] = 10.

  def computational_model( self, s ):
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()

    y = self.solve_model(t_span=[0, t[-1]], y0=y0, N=N, p=p, t_eval=tt)
    S = y[0]
    Idaily = -np.diff(S)

    s['Reference Evaluations'] = list(Idaily)
    s['Dispersion'] = len(Idaily) * [p[-1]]

  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    y = self.solve_model(t_span=[0, t[-1]], y0=y0, N=N, p=p, t_eval=t)
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
    fig = self.new_figure()
    ax  = fig.subplots( 2 )
    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

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


nSamples = 500
x = argparse.Namespace()
x.dataFolder = "data/"
x.country = "switzerland"
x.nPropagation = 100
x.percentages = [0.5]
x.nThreads = 4

a = Model( **vars(x) )

a.sample( nSamples )

a.propagate()

a.save()

a.plot_intervals()
