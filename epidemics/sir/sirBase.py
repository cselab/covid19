#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import requests
import pandas as pd
import io
import os
import numpy as np

from ..epidemics import epidemicsBase
from ..tools.tools import save_file


class sirBase( epidemicsBase ):


  def __init__( self, fileName=[], defaultProperties={}, **kwargs ):

    defaultProperties = { **defaultProperties,
        'country': 'switzerland',
        'populationSize': 8000000,
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
      url = 'https://hgis.uw.edu/virus/assets/virus.csv'
      print(f'Retrieve population data for {self.country} from: {url}')

      s = requests.get(url).content
      df = pd.read_csv(io.StringIO(s.decode('utf-8')))
      if( not self.country in list( df.columns.values ) ):
        sys.exit('Country not in database.')
      d = df[[self.country]].dropna().values.tolist()
      I = [ float(l.split('-')[0]) for  k in d for l in k ]

    N  = len(I)

    self.data['Raw']['Population Size'] = self.populationSize
    self.data['Raw']['Time'] = np.asarray( [ i for i in range(N) ] )
    self.data['Raw']['Infected'] = np.asarray(I)
    self.data['Raw']['Country'] = self.country




  def set_variables_and_distributions( self ):

    self.e['Variables'][0]['Name'] = 'beta'
    self.e['Variables'][0]['Prior Distribution'] = 'Uniform 0'

    self.e['Variables'][1]['Name'] = 'gamma'
    self.e['Variables'][1]['Prior Distribution'] = 'Uniform 1'

    self.e['Variables'][2]['Name'] = '[Sigma]'
    self.e['Variables'][2]['Prior Distribution'] = 'Uniform 2'

    self.nParameters = 3;

    self.e['Distributions'][0]['Name'] = 'Uniform 0'
    self.e['Distributions'][0]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][0]['Minimum'] = 0
    self.e['Distributions'][0]['Maximum'] = +6.0

    self.e['Distributions'][1]['Name'] = 'Uniform 1'
    self.e['Distributions'][1]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][1]['Minimum'] = 0.0
    self.e['Distributions'][1]['Maximum'] = +6.0

    self.e['Distributions'][2]['Name'] = 'Uniform 2'
    self.e['Distributions'][2]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][2]['Minimum'] = 0.0
    self.e['Distributions'][2]['Maximum'] = +600.0




  def sir_rhs( self, t, y, N, beta, gamma ):
    S, I, R = y
    c1 = beta * S * I / N
    c2 = gamma * I
    dSdt = -c1
    dIdt =  c1 - c2
    dRdt =  c2
    return dSdt, dIdt, dRdt




  def solve_ode( y0, T, N, p ):
    sol = solve_ivp( sir_rhs, t_span=[0, T], y0=y0, args=(N, p[0], p[1]), dense_output=True)
    return sol
