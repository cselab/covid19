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



class ModelBase( EpidemicsBase ):


  def __init__( self, **kwargs ):

    self.country        = kwargs.pop('country', 'switzerland')

    self.stdModel       = kwargs.pop('stdModel', 0)
    self.futureDays     = kwargs.pop('futureDays', 2)
    self.nPropagation   = kwargs.pop('nPropagation', 100)
    self.logPlot        = kwargs.pop('logPlot', False)
    self.nValidation    = kwargs.pop('nValidation', 0)
    self.percentages    = kwargs.pop('percentages', [0.5, 0.95, 0.99])

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.country )
    self.propagationData = {}




  # p = [ beta, mu, alpha, Z, D ]
  def seiir_rhs( self, t, y, N, p ):
    S, E, Ir, Iu = y

    c1 = p[0] * S * Ir / N
    c2 = p[1] * p[0] * S * Iu / N
    c3 = p[2] / p[3] * E
    c4 = (1.-p[2]) / p[3] * E


    dSdt  = - c1 - c2
    dEdt  =   c1 + c2 - E/p[3]
    dIrdt =   c3 - Ir/p[4]
    dIudt =   c4 - Iu/p[4]
    return dSdt, dEdt, dIrdt, dIudt




  def solve_ode( self, y0, T, N, p ):
    sol = solve_ivp( self.seiir_rhs, t_span=[0, T], y0=y0, args=(N, p), dense_output=True)
    return sol
