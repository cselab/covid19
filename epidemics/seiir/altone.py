#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

import matplotlib.pyplot as plt
plt.ioff()

from scipy.integrate import solve_ivp, quad

from epidemics.std_models.std_models import *
from epidemics.tools.tools import prepare_folder
from .model_base import *


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'seiir_altone'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data'

    super().__init__( **kwargs )

    self.process_data()




  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )




  def process_data( self ):

    y = self.data['Raw']['Infected']
    t = self.data['Raw']['Time']
    N = self.data['Raw']['Population Size']

    Ir0 = y[0]
    S0  = N - Ir0
    E0  = 0
    Iu0 = 0
    y0  = S0, E0, Ir0, Iu0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = np.diff( y[0:])
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = np.diff( y[0:-self.nValidation] )
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = np.diff( y[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.populationSize
    self.data['Model']['Standard Deviation Model'] = self.stdModel

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1)).tolist()

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )




  def incidence( self, sol, p, t1, t2 ):
    def f(t): y = sol.sol(t); return p[2]/p[3]*y[1];
    return quad(f,t1,t2)[0]




  def computational_model( self, s ):
    p  = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), dense_output=True )
    y = [ self.incidence(sol,p,s-1,s) for s in t ]

    s['Reference Evaluations'] = y
    d = self.data['Model']['y-data']
    s['Standard Deviation Model'] = standard_deviation_models.get( self.stdModel, standardDeviationModelConst)(p[-1],t,d);




  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), dense_output=True )
    y = [ self.incidence(sol,p,s-1,s) for s in t ]

    js = {}
    js['Variables'] = [{}]

    js['Variables'][0]['Name']   = 'Daily Reported Incidence'
    js['Variables'][0]['Values'] = y

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(js['Variables'][0]['Values'])

    d = self.data['Model']['y-data']
    js['Standard Deviation Model'] = standard_deviation_models.get( self.stdModel, standardDeviationModelConst)(p[-1],t,d);

    s['Saved Results'] = js


  def set_variables_for_interval( self ):

    self.intervalVariables = {}

    self.intervalVariables['Daily Reported Incidence'] = {}
    self.intervalVariables['Daily Reported Incidence']['Formula'] = lambda v: v['Daily Reported Incidence']





  def plot_intervals( self ):

    fig = plt.figure(figsize=(12, 8))

    fig.suptitle(self.modelDescription + '  (' + self.country + ')')

    ax  = fig.subplots( 2 )


    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Total Reported Infected(data)', color='black')

    Ns = self.propagatedVariables['Daily Reported Incidence'].shape[0]
    Nt = self.propagatedVariables['Daily Reported Incidence'].shape[1]
    ns = 10;
    samples = np.zeros((Ns*ns,Nt))
    for k in range(Nt):
      x = [ np.random.normal( self.propagatedVariables['Daily Reported Incidence'][:,k],self.propagatedVariables['Standard Deviation'][:,k]) for _ in range(ns) ]
      samples[:,k] = np.asarray(x).flatten()
      samples[:,k] = np.maximum(samples[:,k],0)

    samples = np.cumsum(samples,axis=1)


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

      ax[1].fill_between( self.data['Propagation']['x-data'], q1 , q2,  alpha=0.5, label=f' {100*p:.1f}% credible interval' )


    ax[1].plot( self.data['Propagation']['x-data'], mean, '-', lw=2, label='Mean', color='black')
    ax[1].plot( self.data['Propagation']['x-data'], median, '-', lw=2, label='Median', color='black')


    ax[1].legend(loc='upper left')
    ax[1].set_ylabel( 'Total number of infected' )
    ax[1].set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax[1].grid()


    #----------------------------------------------------------------------------------------------------------------------------------
    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Total Infected(data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

    y = 'Daily Reported Incidence'
    ax[0].plot( self.data['Propagation']['x-data'], self.credibleIntervals[y]['Mean'],   '-', lw=2, label='Mean', color='blue' )
    ax[0].plot( self.data['Propagation']['x-data'], self.credibleIntervals[y]['Median'], '-', lw=2, label='Median', color='black')

    self.credibleIntervals[y]['Intervals'].sort(key = lambda x: x['Percentage'], reverse = True)

    for x in self.credibleIntervals[y]['Intervals']:
      p1 = [ max(k,0) for k in x['Low Interval'] ]
      p2 = x['High Interval']
      p  = 100.*x['Percentage']
      ax[0].fill_between( self.data['Propagation']['x-data'], p1 , p2,  alpha=0.5, label=f' {p:.1f}% credible interval' )

    ax[0].legend(loc='upper left')
    ax[0].set_ylabel( y )
    ax[0].set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax[0].grid()
    if( self.logPlot ): ax[k].set_yscale('log')

    ax[-1].set_xlabel('time in days')

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
