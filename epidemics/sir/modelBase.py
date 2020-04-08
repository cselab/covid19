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


class modelBase( epidemicsBase ):


  def __init__( self, fileName=[], defaultProperties={}, **kwargs ):

    defaultProperties = { **defaultProperties,
        'country': 'switzerland',
        'populationSize': -1,
        'rawData': []
    }

    super().__init__( fileName=fileName, defaultProperties=defaultProperties, **kwargs )

    if fileName == []:
      self.download_raw_data()
      self.propagationData={}




  def download_raw_data( self ):

    if( self.rawData ):
      I = self.rawData
    else:

      if os.path.isfile(self.saveInfo['database']):
          s = load_file(self.saveInfo['database'], 'Downloaded Database', 'pickle')
      else:
        url = 'https://hgis.uw.edu/virus/assets/virus.csv'
        print(f'[Epidemics] Retrieve population data for {self.country} from: {url}')
        s = requests.get(url).content
        save_file( s, self.saveInfo['database'], 'Downloaded Database', 'pickle' )

      df = pd.read_csv(io.StringIO(s.decode('utf-8')))
      if( not self.country in list( df.columns.values ) ):
        sys.exit('Country not in database.')
      d = df[[self.country]].dropna().values.tolist()
      I = [ float(l.split('-')[0]) for  k in d for l in k ]

    N  = len(I)
    if self.populationSize < 0:
      self.populationSize = population_of( self.country )
    self.data['Raw']['Population Size'] = self.populationSize
    self.data['Raw']['Time'] = np.asarray( [ i for i in range(N) ] )
    self.data['Raw']['Infected'] = np.asarray(I)
    self.data['Raw']['Country'] = self.country




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
