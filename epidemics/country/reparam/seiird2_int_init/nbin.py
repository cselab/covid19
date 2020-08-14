import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiird2_int_init.nbin'
    self.modelDescription = 'Fit SEIIRD-2 with Intervention on Daily Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 12
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']), 
            ('Z', *self.defaults['Z']), 
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']),
            ('eps', *self.defaults['eps']),
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            ('e0', *self.defaults['e0']),
            ('iu0', *self.defaults['iu0']),
            ('r', *self.defaults['r'])
            )
    
    return js
