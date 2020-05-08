#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

import matplotlib.pyplot as plt
plt.ioff()

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.tools.tools import prepare_folder, save_file
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'sir_basic_tnrm'
    self.modelDescription = 'Fit SIR on Cummulative Infected Data with Positive Normal likelihood'
    self.likelihoodModel  = 'Positive Normal'

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
    self.e['Distributions'][k]['Minimum'] = 0.00001
    self.e['Distributions'][k]['Maximum'] = 20.




  def computational_model( self, s ):
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.sir_rhs, t_span=[0, t[-1]], y0=y0, args=(N,p), t_eval=t )
    y = ( N - sol.y[0] )

    s['Reference Evaluations'] = y.tolist()
    d = self.data['Model']['y-data']
    s['Standard Deviation'] = ( p[-1] * np.maximum(np.abs(y),1e-4) ).tolist()




  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.sir_rhs, t_span=[0, t[-1]], y0=y0, args=(N,p), t_eval=t )
    y = ( N - sol.y[0] )

    js = {}
    js['Variables'] = [{},{}]

    js['Variables'][0]['Name'] = 'Cummulative Infected'
    js['Variables'][0]['Values'] = y.tolist()

    js['Variables'][1]['Name'] = 'I'
    js['Variables'][1]['Values'] = sol.y[1].tolist()

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = sol.y.shape[1]

    d = self.data['Model']['y-data']
    js['Standard Deviation'] = ( p[-1] * np.maximum(np.abs(y),1e-4) ).tolist()

    s['Saved Results'] = js




  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    ax  = fig.subplots( 1 )

    z = self.data['Model']['y-data']
    ax.plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Total Infected(data)', color='black')

    self.compute_plot_intervals( 'Cummulative Infected', ns, ax, 'Daily Incidence' )

    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)
    plt.show()

    plt.close(fig)
