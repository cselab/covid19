import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.saphire_int.poi'
    self.modelDescription = 'Fit SAPHIRE with Intervention on Daily Infected Data with Poisson likelihood'
    self.likelihoodModel  = 'Poisson'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 9
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
            )
    
    return js
