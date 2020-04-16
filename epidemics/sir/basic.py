#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

import matplotlib.pyplot as plt
plt.ioff()

from scipy.integrate import solve_ivp

from epidemics.std_models.std_models import *
from epidemics.tools.tools import prepare_folder, save_file
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'sir_basic'
    self.modelDescription = 'Fit SIR on Cummulative Infected Data'

    super().__init__( **kwargs )

    self.process_data()


  def process_data( self ):

    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize

    I0 = y[0]
    S0 = N - I0
    y0 = S0, I0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = y[1:]
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = y[1:-self.nValidation]
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = y[-self.nValidation:]

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = N
    self.data['Model']['Standard Deviation Model'] = self.stdModel

    T = t[-1] + self.futureDays
    self.data['Propagation']['x-data'] = np.linspace(0,T,self.nPropagation).tolist()

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )



  def set_variables_and_distributions( self ):

    p = ['beta','gamma','[Sigma]']
    for k,x in enumerate(p):
      self.e['Variables'][k]['Name'] = x
      self.e['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

    self.nParameters = len(p)

    k=0
    self.e['Distributions'][k]['Name'] = 'Prior for beta'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 1.
    self.e['Distributions'][k]['Maximum'] = 50.
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for gamma'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 1.
    self.e['Distributions'][k]['Maximum'] = 50.
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for [Sigma]'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 150.
    self.e['Distributions'][k]['Maximum'] = 800.




  def computational_model( self, s ):
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.sir_rhs, t_span=[0, t[-1]], y0=y0, args=(N,p), t_eval=t )
    y = ( N - sol.y[0]).tolist()


    s['Reference Evaluations'] = y
    d = self.data['Model']['y-data']
    s['Standard Deviation Model'] = standard_deviation_models.get( self.stdModel, standardDeviationModelConst)(p[-1],t,d);




  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.sir_rhs, t_span=[0, t[-1]], y0=y0, args=(N,p), t_eval=t )

    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][0]['Name'] = 'S'
    js['Variables'][0]['Values'] = sol.y[0].tolist()

    js['Variables'].append({})
    js['Variables'][1]['Name'] = 'I'
    js['Variables'][1]['Values'] = sol.y[1].tolist()

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = sol.y.shape[1]

    d = self.data['Model']['y-data']
    js['Standard Deviation Model'] = standard_deviation_models.get( self.stdModel, standardDeviationModelConst)(p[-1],t,d);

    s['Saved Results'] = js




  def set_variables_for_interval( self ):

    self.intervalVariables = {}

    self.intervalVariables['Total Infected'] = {}
    self.intervalVariables['Total Infected']['Formula'] = lambda v: self.regionalData.populationSize - v['S']





  def plot_intervals( self ):

    fig = plt.figure(figsize=(12, 8))

    fig.suptitle(self.modelDescription)

    ax  = fig.subplots( 1 )

    ax.plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Total Infected (data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Total Infected (validation data)', color='black')

    y = 'Total Infected'

    ax.plot( self.data['Propagation']['x-data'], self.credibleIntervals[y]['Mean'],   '-', lw=2, label='Mean', color='blue' )
    ax.plot( self.data['Propagation']['x-data'], self.credibleIntervals[y]['Median'], '-', lw=2, label='Median', color='black')

    self.credibleIntervals[y]['Intervals'].sort(key = lambda x: x['Percentage'], reverse = True)

    for x in self.credibleIntervals[y]['Intervals']:
      p1 = [ max(k,0) for k in x['Low Interval'] ]
      p2 = x['High Interval']
      p  = 100.*x['Percentage']
      ax.fill_between( self.data['Propagation']['x-data'], p1 , p2,  alpha=0.5, label=f' {p:.1f}% credible interval' )

    ax.legend(loc='upper left')
    ax.set_ylabel( y )
    ax.set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax.grid()
    if( self.logPlot ): ax[k].set_yscale('log')

    ax.set_xlabel('time in days')

    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    plt.show()

    plt.close(fig)


    fig = plt.figure(figsize=(12, 8))
    ax  = fig.subplots(1)

    R0 = self.parameters[0]['Values'] / self.parameters[1]['Values']

    ax.hist( R0 , 100, density=True, facecolor='g', alpha=0.75)
    ax.set_xlabel('R0')
    ax.grid()

    file = os.path.join(self.saveInfo['figures'],'R0.png')
    prepare_folder( os.path.dirname(file), clean=False )
    fig.savefig(file)
