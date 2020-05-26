#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   31/3/2020
# Email:  garampat@ethz.ch

import requests
import pandas as pd
import io
import os
import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
plt.ioff()

from epidemics.tools.tools import prepare_folder, save_file
from epidemics.data.combined import RegionalData
from epidemics.epidemics import EpidemicsBase



class ModelBase( EpidemicsBase ):


  def __init__( self, **kwargs ):

    self.country        = kwargs.pop('country', 'switzerland')

    self.stdModel       = kwargs.pop('stdModel', 0)
    self.futureDays     = kwargs.pop('futureDays', 2)
    self.nPropagation   = kwargs.pop('nPropagation', 100)
    self.logPlot        = kwargs.pop('logPlot', False)
    self.nValidation    = kwargs.pop('nValidation', 0)
    self.percentages    = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess     = kwargs.pop('preprocess', False)

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.country )
    self.propagationData = {}





  def process_data( self ):

    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize

    i0 = 0

    if i0==0:
      Ir0 = y[0]
    else:
      Ir0 = y[i0]-y[i0-1]

    S0  = N - Ir0
    E0  = 0
    Iu0 = 0
    y0  = S0, E0, Ir0, Iu0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[i0+1:]
      self.data['Model']['y-data'] = np.diff(y[i0:])
    else:
      self.data['Model']['x-data'] = t[i0+1:-self.nValidation]
      self.data['Model']['y-data'] = np.diff( y[i0:-self.nValidation] )
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = np.diff( y[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize
    self.data['Model']['Standard Deviation Model'] = self.stdModel

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(i0,T,int(T+1))

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )




  #       0    1     2   3  4    5     6
  # p = [ R0, mu, alpha, Z, D, delta, td ]
  def seiir_rhs( self, t, y, N, p ):
    S, E, Ir, Iu = y

    beta = p[0]*p[4]

    if t>p[6] :
      beta = p[0]*p[5]

    c1 = beta * S * Ir / N
    c2 = p[1] * beta * S * Iu / N
    c3 = p[2] / p[3] * E
    c4 = (1.-p[2]) / p[3] * E


    dSdt  = - c1 - c2
    dEdt  =   c1 + c2 - E/p[3]
    dIrdt =   c3 - Ir/p[4]
    dIudt =   c4 - Iu/p[4]
    return dSdt, dEdt, dIrdt, dIudt




  def solve_ode( self, y0, T, N, p ):
    sol = solve_ivp( self.seiir_rhs, t_span=[0, T], y0=y0, args=(N, p), dense_output=True)
    return sol




  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    ax  = fig.subplots( 2 )

    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected (data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

    self.compute_plot_intervals( 'Daily Reported Incidence', ns, ax[0], 'Daily Reported Incidence' )

    #----------------------------------------------------------------------------------------------------------------------------------
    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected(data)', color='black')

    self.compute_plot_intervals( 'Daily Reported Incidence', ns, ax[1], 'Cummulative number of reported infected', cummulate=1)

    #----------------------------------------------------------------------------------------------------------------------------------

    ax[-1].set_xlabel('time in days')

    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    plt.show()

    plt.close(fig)

    #----------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------

    fig = self.new_figure()

    ax  = fig.subplots( 1 )

    if self.parameters[0]['Name'] == 'R0':
      ax.hist( self.parameters[0]['Values'], bins=40, density=1)
    else:
      ax.hist( self.parameters[0]['Values']/self.parameters[1]['Values'], bins=40, density=1)

    file = os.path.join(self.saveInfo['figures'],'R0.png');
    fig.savefig(file)

    plt.show()

    plt.close(fig)
