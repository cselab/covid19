#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

import matplotlib.pyplot as plt
plt.ioff()

from scipy.integrate import solve_ivp, quad

from  .modelBase import *
from ..std_models.std_models import *
from ..tools.tools import prepare_folder


class model( modelBase ):


  def __init__( self, fileName=[], **kwargs ):

    self.modelName        = 'seiird_altone'
    self.modelDescription = 'Fit SEIIRD on Daily Infected and Deaths Data'

    defaultProperties = {
        'stdModel': 0,
        'futureDays': 2,
        'nPropagation': 100,
        'logPlot': False,
        'nValidation': 0
    }

    super().__init__( fileName=fileName, defaultProperties=defaultProperties, **kwargs )

    if fileName == []:
      self.process_data()




  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )




  def process_data( self ):

    i = self.data['Raw']['Infected']
    d = self.data['Raw']['Dead']
    t = self.data['Raw']['Time']
    N = self.data['Raw']['Population Size']

    Ir0 = i[0]
    D0  = d[0]
    S0  = N - Ir0 - D0
    E0  = 0
    Iu0 = 0
    y0  = S0, E0, Ir0, Iu0, D0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      I = np.diff( i[0:])
      D = np.diff( d[0:])
      self.data['Model']['I-data'] = I
      self.data['Model']['D-data'] = D
      self.data['Model']['y-data'] = np.concatenate((I,D))
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      I = np.diff( i[0:-self.nValidation] )
      D = np.diff( d[0:-self.nValidation] )
      self.data['Model']['I-data'] = I
      self.data['Model']['D-data'] = D
      self.data['Model']['y-data'] = np.concatenate((I,D))

      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['I-data'] = np.diff( i[-self.nValidation-1:] )
      self.data['Validation']['D-data'] = np.diff( d[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.populationSize
    self.data['Model']['Standard Deviation Model'] = self.stdModel

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1)).tolist()

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )



  def infected_incidence( self, sol, p, t1, t2 ):
    def f(t): y = sol.sol(t); return p[2]/p[3]*y[1];
    return quad(f,t1,t2)[0]




  def deaths_incidence( self, sol, p, t1, t2 ):
    def f(t): y = sol.sol(t); return p[5]*y[2];
    return quad(f,t1,t2)[0]




  def computational_model( self, s ):
    p  = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), dense_output=True )
    I = [ self.infected_incidence(sol,p,s-1,s) for s in t ]
    D = [ self.deaths_incidence(sol,p,s-1,s) for s in t ]

    s['Reference Evaluations'] = I + D
    d = self.data['Model']['y-data']
    std_model = standard_deviation_models.get( self.stdModel, standardDeviationModelError)(p[-1],t,d);
    s['Standard Deviation Model'] = 2*std_model




  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), dense_output=True )
    I = [ self.infected_incidence(sol,p,s-1,s) for s in t ]
    D = [ self.deaths_incidence(sol,p,s-1,s) for s in t ]

    js = {}
    js['Variables'] = [{},{},{}]

    js['Variables'][0]['Name']   = 'Daily Reported Infected Incidence'
    js['Variables'][0]['Values'] = I

    js['Variables'][1]['Name']   = 'Daily Deaths Incidence'
    js['Variables'][1]['Values'] = D

    js['Variables'][2]['Name']   = 'Daily Infected'
    js['Variables'][2]['Values'] = D

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(js['Variables'][0]['Values'])

    d = self.data['Model']['y-data']
    js['Standard Deviation Model'] = standard_deviation_models.get( self.stdModel, standardDeviationModelError)(p[-1],t,d);

    s['Saved Results'] = js




  def set_variables_for_interval( self ):

    self.intervalVariables = {}

    self.intervalVariables['Daily Reported Infected Incidence'] = {}
    self.intervalVariables['Daily Reported Infected Incidence']['Formula'] = lambda v: v['Daily Reported Infected Incidence']

    self.intervalVariables['Daily Deaths Incidence'] = {}
    self.intervalVariables['Daily Deaths Incidence']['Formula'] = lambda v: v['Daily Deaths Incidence']

    self.intervalVariables['Daily Infected'] = {}
    self.intervalVariables['Daily Infected']['Formula'] = lambda v: v['Daily Infected']




  def plot_variable(self, ax, varName ):

    ax.plot( self.data['Propagation']['x-data'], self.credibleIntervals[varName]['Mean'],   '-', lw=2, label='Mean', color='blue' )
    ax.plot( self.data['Propagation']['x-data'], self.credibleIntervals[varName]['Median'], '-', lw=2, label='Median', color='black')

    self.credibleIntervals[varName]['Intervals'].sort(key = lambda x: x['Percentage'], reverse = True)

    for x in self.credibleIntervals[varName]['Intervals']:
      p1 = [ max(k,0) for k in x['Low Interval'] ]
      p2 = x['High Interval']
      p  = 100.*x['Percentage']
      ax.fill_between( self.data['Propagation']['x-data'], p1 , p2,  alpha=0.5, label=f' {p:.1f}% credible interval' )

    ax.legend(loc='upper left')
    ax.set_ylabel( varName )
    ax.set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax.grid()
    if( self.logPlot ): ax[k].set_yscale('log')




  def plot_cummulative_variable( self, ax, varName, label ):

    Ns = self.propagatedVariables[varName].shape[0]
    Nt = self.propagatedVariables[varName].shape[1]
    ns = 10;
    samples = np.zeros((Ns*ns,Nt))
    for k in range(Nt):
      x = [ np.random.normal( self.propagatedVariables[varName][:,k],self.propagatedVariables['Standard Deviation'][:,k]) for _ in range(ns) ]
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

      ax.fill_between( self.data['Propagation']['x-data'], q1 , q2,  alpha=0.5, label=f' {100*p:.1f}% credible interval' )


    ax.plot( self.data['Propagation']['x-data'], mean, '-', lw=2, label='Mean', color='black')
    ax.plot( self.data['Propagation']['x-data'], median, '-', lw=2, label='Median', color='black')


    ax.legend(loc='upper left')
    ax.set_ylabel( label )
    ax.set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
    ax.grid()





  def plot_intervals( self ):

    fig = plt.figure(figsize=(12, 8))

    fig.suptitle(self.modelDescription + '  (' + self.country + ')')

    ax = fig.subplots(2,2)

    #----------------------------------------------------------------------------------------------------------------------------------
    # Daily plots
    ax[0,0].plot( self.data['Model']['x-data'], self.data['Model']['I-data'], 'o', lw=2, label='Daily Reported Infected (data)', color='black')

    if self.nValidation > 0:
      ax[0,0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='(validation data)', color='black')

    self.plot_variable(ax[0,0], 'Daily Reported Infected Incidence' )


    ax[1,0].plot( self.data['Model']['x-data'], self.data['Model']['D-data'], 'o', lw=2, label='Daily Deaths (data)', color='black')

    if self.nValidation > 0:
      ax[1,0].plot( self.data['Validation']['x-data'], self.data['Validation']['D-data'], 'x', lw=2, label='(validation data)', color='black')

    self.plot_variable(ax[1,0], 'Daily Deaths Incidence' )

    ax[-1,0].set_xlabel('time in days')


    #----------------------------------------------------------------------------------------------------------------------------------
    # Cummulative plots
    z = np.cumsum(self.data['Model']['I-data'])
    ax[0,1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Total Reported Infected (data)', color='black')

    self.plot_cummulative_variable( ax[0,1], 'Daily Reported Infected Incidence', 'Total Infected' )

    z = np.cumsum(self.data['Model']['D-data'])
    ax[1,1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Total Deaths (data)', color='black')

    self.plot_cummulative_variable( ax[1,1], 'Daily Deaths Incidence', 'Total Deaths' )

    ax[-1,1].set_xlabel('time in days')

    #----------------------------------------------------------------------------------------------------------------------------------
    # Save figure
    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    plt.show()

    plt.close(fig)
