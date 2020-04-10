#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import requests
import pandas as pd
import io
import os
import numpy as np
from scipy.integrate import solve_ivp

from ..epidemics import epidemicsBase
from ..tools.tools import save_file, load_file
from ..tools.population_of import population_of
from ..tools.database import regionalData


class modelBase( epidemicsBase ):


  def __init__( self, fileName=None, defaultProperties={}, **kwargs ):

    defaultProperties = { **defaultProperties,
        'country': 'switzerland'
    }

    super().__init__( fileName=fileName, defaultProperties=defaultProperties, **kwargs )

    if not fileName:
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
