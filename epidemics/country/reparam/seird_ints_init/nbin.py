import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seird_ints_init.nbin'
    self.modelDescription = 'Fit SEIRD with Intervention on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 8
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Z', *self.defaults['Z']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('e0', *self.defaults['e0']),
            ('r', *self.defaults['r']),
            )
    
    return js
