#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import requests
import io
import os
import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
plt.ioff()

from epidemics.data.combined import RegionalData
from epidemics.epidemics import EpidemicsBase
import epidemics.ode_solver as solver
from epidemics.utils.misc import prepare_folder, save_file

class ModelBase( EpidemicsBase ):


  def __init__( self, **kwargs ):

    self.country      = kwargs.pop('country', 'switzerland')
    self.futureDays   = kwargs.pop('futureDays', 2)
    self.nPropagation = kwargs.pop('nPropagation', 100)
    self.logPlot      = kwargs.pop('logPlot', False)
    self.nValidation  = kwargs.pop('nValidation', 0)
    self.percentages  = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess   = kwargs.pop('preprocess', False)

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.country,self.preprocess)
    self.propagationData = {}




  def save_data_path( self ):
      if hasattr(self, 'property'):
        return ( self.dataFolder, self.country, self.modelName )
      else:
        return ( self.dataFolder, self.country )



  # R0, gamma
  def sir_rhs( self, t, y, N, p ):
      S, I, S0, I0, S1, I1 = y

      SI = S*I
      p01N = p[0] * p[1] / N

      c1 = p01N*SI
      c2 = p[1]*I
      c3 = p[1]/N*SI
      c4 = p[0]/N*SI

      a11 = -p01N*I; a12 = -p01N*S
      a21 =  p01N*I; a22 =  p01N*S - p[1]

      dSdt = -c1
      dIdt =  c1 - c2

      dS0dt = a11*S0 + a12*I0 - c3
      dI0dt = a21*S0 + a22*I0 + c3

      dS1dt = a11*S1 + a12*I1 - c4
      dI1dt = a21*S1 + a22*I1 + c4 - I

      return dSdt, dIdt, dS0dt, dI0dt, dS1dt, dI1dt




  def solve_ode( self, y0, T, N, p ):
    sol = solve_ivp( self.sir_rhs, t_span=[0, T], y0=y0, args=(N, p), dense_output=True)
    return sol




  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    ax  = fig.subplots( 2 )

    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', ns, ax[0], 'Daily Incidence' )

    #----------------------------------------------------------------------------------------------------------------------------------
    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected(data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', ns, ax[1], 'Cummulative number of infected', cummulate=1)

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
