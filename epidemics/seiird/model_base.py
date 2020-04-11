#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   31/3/2020
# Email:  garampat@ethz.ch

import requests
import pandas as pd
import io
import os
import numpy as np
from scipy.integrate import solve_ivp


from epidemics.data.combined import RegionalData
from epidemics.epidemics import EpidemicsBase
from epidemics.tools.tools import save_file, load_file


class ModelBase( EpidemicsBase ):


  def __init__( self, **kwargs ):

    self.country = kwargs.pop('country', 'switzerland')

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.saveInfo['database'], self.country )
    self.propagationData={}




  def set_variables_and_distributions( self ):

    p = [ 'beta', 'mu', 'alpha', 'Z', 'D', 'd', '[Sigma]' ]

    for k,x in enumerate(p):
      self.e['Variables'][k]['Name'] = x
      self.e['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

    self.nParameters = len(p)

    k=0
    self.e['Distributions'][k]['Name'] = 'Prior for beta'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0
    self.e['Distributions'][k]['Maximum'] = 5
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for mu'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0
    self.e['Distributions'][k]['Maximum'] = 1
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for alpha'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0.1
    self.e['Distributions'][k]['Maximum'] = 1.0
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for Z'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0.1
    self.e['Distributions'][k]['Maximum'] = 5
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for D'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 2
    self.e['Distributions'][k]['Maximum'] = 5
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for d'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 0.1
    self.e['Distributions'][k]['Maximum'] = 10
    k+=1

    self.e['Distributions'][k]['Name'] = 'Prior for [Sigma]'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 10
    self.e['Distributions'][k]['Maximum'] = 1000




  # p = [ beta, mu, alpha, Z, D, d ]
  def seiir_rhs( self, t, y, N, p ):
    S, E, Ir, Iu, D = y

    c1 = p[0] * S * Ir / N
    c2 = p[1] * p[0] * S * Iu / N
    c3 = p[2] / p[3] * E
    c4 = (1.-p[2]) / p[3] * E
    c5 = p[5] / p[4] * Ir

    dSdt  = - c1 - c2
    dEdt  =   c1 + c2 - (c3 + c4)
    dIrdt =   c3 - Ir/p[4] - c5
    dIudt =   c4 - Iu/p[4]
    dDdt  =   c5
    return dSdt, dEdt, dIrdt, dIudt, dDdt




  def solve_ode( self, y0, T, N, p ):
    sol = solve_ivp( self.seiir_rhs, t_span=[0, T], y0=y0, args=(N, p), dense_output=True)
    return sol
