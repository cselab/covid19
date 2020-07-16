#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

import matplotlib.pyplot as plt
plt.ioff()

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.utils.misc import prepare_folder, save_file
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'seiir.tnrm'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data with Positive Normal likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

    self.process_data()




  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )




  def get_variables_and_distributions( self ):
    p = [ 'R0', 'mu', 'alpha', 'Z', 'D', '[Sigma]' ]

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
    js['Distributions'][k]['Minimum'] = 0.1
    js['Distributions'][k]['Maximum'] = 2.5

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for mu'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0
    js['Distributions'][k]['Maximum'] = 0.1

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for alpha'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.
    js['Distributions'][k]['Maximum'] = 0.1

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for Z'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0
    js['Distributions'][k]['Maximum'] = 40

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for D'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 2000
    js['Distributions'][k]['Maximum'] = 4000

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for [Sigma]'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.01
    js['Distributions'][k]['Maximum'] = 20

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




  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solve_ivp( self.seiir_rhs, t_span=[0, t[-1]], y0=y0, args=(N, p), t_eval=t )

    y   = - p[2] * ( np.diff(sol.y[0]) + np.diff(sol.y[1]) )
    std = ( p[-1] * np.maximum(np.abs(y),1e-4) )

    y   = [float(y0[2])] + y.tolist()
    std = [1e-5] + std.tolist()

    js = {}
    js['Variables'] = [{}]

    js['Variables'][0]['Name'] = 'Daily Reported Incidence'
    js['Variables'][0]['Values'] = y

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(js['Variables'][0]['Values'])

    js['Standard Deviation'] = std

    s['Saved Results'] = js
