import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiir_int_nogamma.nbin'
    self.modelDescription = 'Fit SEIIR with Interventions on Daily Infected Data with Negative Binomial Likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )

  def get_variables_and_distributions( self ):
 
    self.nParameters = 8
    js = self.get_uniform_priors(
            ('R0', 1.0, 10.0), 
            ('Z', 1.0, 30.), 
            ('mu', 0.0, 1.0), 
            ('alpha', 0., 1.0),
            ('tact', 0.0, 100.),
            ('dtact', 9.99, 10.01),
            ('kbeta', 0.1, 1.0),
            ('r', 1e-6, 1.0),
            )
    
    return js
