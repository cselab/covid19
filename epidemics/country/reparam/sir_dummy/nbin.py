import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir_dummy.nbin'
    self.modelDescription = 'Fit SIR with Intervention on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 5
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('D', *self.defaults['D']), 
            ('tact', *self.defaults['tact']), 
            ('dtact', *self.defaults['dtact']), 
            ('r', *self.defaults['r']),
            )
    
    return js
