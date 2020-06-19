import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir_int.tnrm'
    self.modelDescription = 'Fit SIR with interventions on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

 
  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('R0', 1.0, 10.0), 
            ('D', 1.0, 30.), 
            ('tact', 0, 100),
            ('dtact', 9.99, 10.01),
            ('kbeta', 0.1, 1.0),
            ('Sigma', 1e-6, 100)
            )
    
    return js
