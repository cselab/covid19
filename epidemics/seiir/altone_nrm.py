#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

import matplotlib.pyplot as plt
plt.ioff()

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.std_models.std_models import *
from epidemics.tools.tools import prepare_folder, save_file
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'seiir_altone_nrm'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data with Normal likelihood'
    self.likelihoodModel  = 'Normal'

    super().__init__( **kwargs )

    self.process_data()




  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )




  def get_variables_and_distributions( self ):
    p = [ 'beta', 'mu', 'alpha', 'Z', 'D', '[Sigma]' ]

    js = {}
    js['Variables']=[]
    js['Distributions']=[]

    for k,x in enumerate(p):
      js['Variables'].append({})
      js['Variables'][k]['Name'] = x
      js['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

    self.nParameters = len(p)

    k=0
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for beta'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0
    js['Distributions'][k]['Maximum'] = 100

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for mu'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0
    js['Distributions'][k]['Maximum'] = 1

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for alpha'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.1
    js['Distributions'][k]['Maximum'] = 1.0

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for Z'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0
    js['Distributions'][k]['Maximum'] = 10

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for D'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0
    js['Distributions'][k]['Maximum'] = 100

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for [Sigma]'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.01
    js['Distributions'][k]['Maximum'] = 10

    return js




  def computational_model( self, s ):
    p  = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), t_eval=tt )

    y = - p[2] * ( np.diff(sol.y[0]) + np.diff(sol.y[1]) )

    s['Reference Evaluations'] = y.tolist()
    s['Standard Deviation'] = ( p[-1] * np.maximum(np.abs(y),1e-4) ).tolist()
    # s['Standard Deviation'] = ( p[-1] * np.minimum( np.maximum(np.abs(y),1e-4), 2e3) ).tolist()
    # s['Standard Deviation'] = len(y)*[ 100*p[-1]]


  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), t_eval=t )

    y = - p[2] * ( np.diff(sol.y[0]) + np.diff(sol.y[1]) )

    y = [float(y0[2])] + y.tolist()


    js = {}
    js['Variables'] = [{}]

    js['Variables'][0]['Name'] = 'Daily Reported Incidence'
    js['Variables'][0]['Values'] = y

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(js['Variables'][0]['Values'])

    js['Standard Deviation'] = ( p[-1] * np.maximum(np.abs(y),1e-4) ).tolist()
    # js['Standard Deviation'] = ( p[-1] * np.minimum( np.maximum(np.abs(y),1e-4), 2e3) ).tolist()
    # js['Standard Deviation'] = len(y)*[100*p[-1]]

    s['Saved Results'] = js




  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    ax  = fig.subplots( 2 )

    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

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
