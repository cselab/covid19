import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam_beta.seiir.tnrm'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('beta', 0.0, 10.0),
            ('D', 1.0, 30.), 
            ('Z', 1.0, 30.0), 
            ('mu', 0.0, 1.0), 
            ('alpha', 0., 1.0),
            ('Sigma', 1e-6, 100)
            )
    
    return js


