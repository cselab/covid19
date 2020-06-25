import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):

  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir.nbin'
    self.modelDescription = 'Fit SIR on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):

    self.nParameters = 3
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('D', *self.defaults['D']), 
            ('r', *self.defaults['r']))
    
    return js
