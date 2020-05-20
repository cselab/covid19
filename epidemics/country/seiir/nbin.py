#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.tools.tools import prepare_folder, save_file
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.seiir.nbin'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )

    self.process_data()




  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )




  def get_variables_and_distributions( self ):
    p = [ 'beta', 'mu', 'alpha', 'Z', 'D', '[r]' ]

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
    js['Distributions'][k]['Minimum'] = 0.1
    js['Distributions'][k]['Maximum'] = 100.

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
    js['Distributions'][k]['Maximum'] = 50

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for D'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1000
    js['Distributions'][k]['Maximum'] = 8000

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for [r]'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.001
    js['Distributions'][k]['Maximum'] = 10

    return js


  def computational_model( self, s ):
 
    p = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt,N=N,p=p)

    # get incidents
    y = -np.diff(sol.y[0])
     
    eps = 1e-32
    y[y < eps] = eps
  
    if(self.sampler == 'mTMCMC'):
        print("[Epidemics] mTMCMC not yet available for nbin")
        sys.exit(0)

    s['Reference Evaluations'] = list(y)
    s['Dispersion'] = ( p[-1] * y ).tolist()

  def computational_model_propagate( self, s ):
    p = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']
 
    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1],t_eval=t.tolist(), N=N,p=p)
    
    y = -np.diff(sol.y[0])
    y = np.append(0, y)

    eps = 1e-32
    y[y < eps] = eps
 
    js = {}
    js['Variables'] = [{}]

    js['Variables'][0]['Name'] = 'Daily Reported Incidence'
    js['Variables'][0]['Values'] = list(y)

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    js['Dispersion'] = ( p[-1] * y ).tolist()

    s['Saved Results'] = js
