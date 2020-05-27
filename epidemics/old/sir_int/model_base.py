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
from epidemics.tools.tools import prepare_folder, save_file


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
    return ( self.dataFolder, self.country, self.modelName )



  # beta, gamma, delta, t_a
  def sir_rhs( self, t, y, N, p ):

    factor = 0.5
    t_start = 45
    t_end = 55

    S, I = y

    R0 = p[0]
    gamma = p[1]

    if( t<p[3] ):
        beta = R0*gamma

    elif t>=p[3] and t< t_start:
        beta = R0*gamma*p[2]

    elif t >= t_start and t< t_end:
        frac = (t-t_start)/(t_end-t_start)
        beta = R0*gamma*p[2]*(1+factor*frac)
    else:
        beta = R0*gamma*p[2]*(1+factor )

    c1 = beta * S * I / N
    c2 = gamma * I
    dSdt = -c1
    dIdt =  c1 - c2
    return dSdt, dIdt




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
    ax[0].set_yscale('log')

    #----------------------------------------------------------------------------------------------------------------------------------
    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected (data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', ns, ax[1], 'Cummulative number of infected', cummulate=1)

    #----------------------------------------------------------------------------------------------------------------------------------
    ax[-1].set_xlabel('time in days')
    ax[-1].set_yscale('log')
    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    plt.show()

    plt.close(fig)

    #----------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------

    fig = self.new_figure()

    ax  = fig.subplots( 2 )

    if self.parameters[0]['Name'] == 'R0':
      R0 = self.parameters[0]['Values']
    else:
      R0 = self.parameters[0]['Values']/self.parameters[1]['Values']

    ax[0].hist( R0, 40, density=1)
    ax[1].hist( R0*self.parameters[2]['Values'], 40, density=1)

    file = os.path.join(self.saveInfo['figures'],'R0.png');
    fig.savefig(file)

    plt.show()

    plt.close(fig)
