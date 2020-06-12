import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiir.nrm'
    self.modelDescription = 'Fit SEIIR on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Normal'

    super().__init__( **kwargs )

  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('R0', 0.0, 2.0), 
            ('mu', 0.0, 2.0), 
            ('alpha', 0., 1.0),
            ('Z', 0, 15.0), 
            ('D', 1.0, 30.), 
            ('Sigma', 1e-6, 100)
            )
    
    return js



