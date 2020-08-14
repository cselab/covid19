import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.saphire_int.geo'
    self.modelDescription = 'Fit SAPHIRE with Intervention on Daily Infected Data with Geometric likelihood'
    self.likelihoodModel  = 'Geometric'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 9
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Zl', *self.defaults['Zl']),
            ('Y', *self.defaults['Y']), 
            ('mu', *self.defaults['alpha']), 
            ('alpha', *self.defaults['alpha']), 
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            )
    
    return js
