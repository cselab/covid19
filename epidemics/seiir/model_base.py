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

from epidemics.data.population import get_country_population
from epidemics.epidemics import EpidemicsBase
from epidemics.tools.tools import save_file, load_file


class ModelBase( EpidemicsBase ):


  def __init__( self, defaultProperties={}, **kwargs ):

    defaultProperties = { **defaultProperties,
        'country': 'switzerland',
        'populationSize': -1,
        'rawData': []
    }

    super().__init__( defaultProperties=defaultProperties, **kwargs )

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
      self.populationSize = get_country_population( self.country )
    self.data['Raw']['Population Size'] = self.populationSize
    self.data['Raw']['Time'] = np.asarray( [ i for i in range(N) ] )
    self.data['Raw']['Infected'] = np.asarray(I)
    self.data['Raw']['Country'] = self.country




  def set_variables_and_distributions( self ):

    p = [ 'beta', 'mu', 'alpha', 'Z', 'D', '[Sigma]' ]

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

    self.e['Distributions'][k]['Name'] = 'Prior for [Sigma]'
    self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
    self.e['Distributions'][k]['Minimum'] = 10
    self.e['Distributions'][k]['Maximum'] = 1000




  # p = [ beta, mu, alpha, Z, D ]
  def seiir_rhs( self, t, y, N, p ):
    S, E, Ir, Iu = y

    c1 = p[0] * S * Ir / N
    c2 = p[1] * p[0] * S * Iu / N
    c3 = p[2] / p[3] * E
    c4 = (1.-p[2]) / p[3] * E


    dSdt  = - c1 - c2
    dEdt  =   c1 + c2 - (c3 + c4)
    dIrdt =   c3 - Ir/p[4]
    dIudt =   c4 - Iu/p[4]
    return dSdt, dEdt, dIrdt, dIudt




  def solve_ode( self, y0, T, N, p ):
    sol = solve_ivp( self.seiir_rhs, t_span=[0, T], y0=y0, args=(N, p), dense_output=True)
    return sol
