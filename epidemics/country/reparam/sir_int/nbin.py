import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir_int.nbin'
    self.modelDescription = 'Fit SIR with Intervention on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('R0', 0.0, 2.0), 
            ('D', 1.0, 30.), 
            ('tact', 1, 80),
            ('dtact', 0.0, 30),
            ('kbeta', 0.1, 10),
            ('[r]', 0.01, 100),
            )
    
    return js