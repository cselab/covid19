import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam_beta.sir_intexp.nbin'
    self.modelDescription = 'Fit SIR with Intervention on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('beta', 0.0, 10.0), 
            ('D', 1.0, 30.), 
            ('tact', 0, 100),
            ('k', 0.0, 1.0),
            ('r', 1e-6, 1.0),
            )
    
    return js
