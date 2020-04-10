#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import requests
import io
import os
import numpy as np
from scipy.integrate import solve_ivp

from epidemics.data.combined import RegionalData
from epidemics.epidemics import EpidemicsBase


class ModelBase( EpidemicsBase ):


  def __init__( self, **kwargs ):

    self.country      = kwargs.pop('country', 'switzerland')
    self.stdModel     = kwargs.pop('stdModel', 0)
    self.futureDays   = kwargs.pop('futureDays', 2)
    self.nPropagation = kwargs.pop('nPropagation', 100)
    self.logPlot      = kwargs.pop('logPlot', False)
    self.nValidation  = kwargs.pop('nValidation', 0)

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.country )
    self.propagationData = {}




  def save_data_path( self ):
    return ( self.dataFolder, self.country, self.modelName )




  def sir_rhs( self, t, y, N, beta, gamma ):
    S, I = y
    c1 = beta * S * I / N
    c2 = gamma * I
    dSdt = -c1
    dIdt =  c1 - c2
    return dSdt, dIdt




  def solve_ode( self, y0, T, N, p ):
    sol = solve_ivp( self.sir_rhs, t_span=[0, T], y0=y0, args=(N, p[0], p[1]), dense_output=True)
    return sol
