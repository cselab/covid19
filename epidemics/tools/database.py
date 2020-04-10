import os
import io
import requests
import pandas as pd
import numpy as np

from .tools import load_file, save_file
from .population_of import population_of


class regionalData( ):


  def __init__( self, dbFileName, country  ):

    if os.path.isfile( dbFileName ):
        s = load_file( dbFileName, 'Downloaded Database', 'pickle')
    else:
      url = 'https://hgis.uw.edu/virus/assets/virus.csv'
      print(f'[Epidemics] Fetch epidemics data for {country} from: {url}')
      s = requests.get(url).content
      save_file( s, dbFileName, 'Downloaded Database', 'pickle' )

    df = pd.read_csv(io.StringIO(s.decode('utf-8')))
    if( not country in list( df.columns.values ) ):
      sys.exit('[Epidemics] Country not in database.')
    d = df[[country]].dropna().values.tolist()
    I = [ float(k[0].split('-')[0]) for  k in d ]
    D = [ float(k[0].split('-')[3]) for  k in d ]

    N  = len(I)

    self.country  = country
    self.populationSize = population_of( country )
    self.time     = np.asarray( [ i for i in range(N) ] )
    self.infected = np.asarray(I)
    self.deaths     = np.asarray(D)
