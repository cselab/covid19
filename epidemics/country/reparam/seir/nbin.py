import numpy as np
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seir.nbin'
    self.modelDescription = 'Fit SEIR on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 4
    js = self.get_uniform_priors(
            ('R0', 0.5, 5.0), 
            ('D', 0.0, 30.0), 
            ('Z', 0.0, 30.0), 
            ('r', 1e-6, 1.0),
            )
    
    return js
