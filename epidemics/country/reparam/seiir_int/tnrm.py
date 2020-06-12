import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiir_int.tnrm'
    self.modelDescription = 'Fit SEIIR with Interventions on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

  def get_variables_and_distributions( self ):
 
    self.nParameters = 9
    js = self.get_uniform_priors(
            ('R0', 0.5, 2.0), 
            ('D', 0.0, 30.),  
            ('Z', 0.0, 30.), 
            ('mu', 0.0, 1.0), 
            ('alpha', 0., 1.0),
            ('tact', 0.0, 100.),
            ('dtact', 0.0, 50.),
            ('kbeta', 0.0, 1.0),
            ('Sigma', 1e-6, 100)
            )
    
    return js
