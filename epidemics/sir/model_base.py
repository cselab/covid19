#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import requests
import io
import os
import numpy as np
from scipy.integrate import solve_ivp

from epidemics.epidemics import EpidemicsBase
from epidemics.tools.tools import save_file, load_file
from epidemics.tools.database import regionalData


class ModelBase( EpidemicsBase ):


  def __init__( self, defaultProperties={}, **kwargs ):

    defaultProperties = { **defaultProperties,
        'country': 'switzerland'
    }

    super().__init__( defaultProperties=defaultProperties, **kwargs )

    self.regionalData = regionalData( self.saveInfo['database'], self.country )
    self.propagationData={}




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
