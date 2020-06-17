import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):

  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam_beta.sir.tnrm'
    self.modelDescription = 'Fit SIR on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 3
    js = self.get_uniform_priors(
            ('beta', 0.0, 10.0), 
            ('D', 1.0, 30.0), 
            ('Sigma', 1e-6, 100)
            )
    
    return js
