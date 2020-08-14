import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.saphire_int_init.nbin'
    self.modelDescription = 'Fit SAPHIRE with Intervention on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 13
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Zl', *self.defaults['Zl']),
            ('Y', *self.defaults['Y']), 
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']), 
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            ('e0', *self.defaults['e0']),
            ('p0', *self.defaults['p0']),
            ('iu0', *self.defaults['iu0']),
            ('r', *self.defaults['r']),
            )
    
    return js
