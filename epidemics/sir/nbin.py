#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.tools.tools import save_file
from .model_base import ModelBase
import epidemics.ode_solver as solver

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'sir.nbin'
    self.modelDescription = 'Fit SIR on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )

    self.process_data()




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
      print('Data t {} y {}'.format(len(t[1:]),len(y[0:])))
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




  def get_variables_and_distributions( self ):

    p = ['R0','gamma','[r]']

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
    js['Distributions'][k]['Name'] = 'Prior for R0'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1.
    js['Distributions'][k]['Maximum'] = 30.

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for gamma'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1.
    js['Distributions'][k]['Maximum'] = 30.

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for [r]'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.01
    js['Distributions'][k]['Maximum'] = 10.

    return js




  def computational_model( self, s ):
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = solver.solve_ode(self.sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = tt,backend=self.backend)    
    y = -(sol.y[0][1:]-sol.y[0][:-1])
    
    # Need to take pytorch gradients before here
    y = solver.to_list(y)
    
    s['Reference Evaluations'] = y
    s['Dispersion'] = len(y)*[p[-1]]



  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solver.solve_ode(self.sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = t.tolist(),backend='numpy')
    y = -(sol.y[0][1:]-sol.y[0][:-1])
    y = solver.append_zero(y)
    
    # Need to take pytorch gradients before here
    y = solver.to_list(y)

    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][0]['Name']   = 'Daily Incidence'
    js['Variables'][0]['Values'] = y

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    js['Dispersion'] = (len(y))*[p[-1]]

    s['Saved Results'] = js
