import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiir_ints.tstudent_alt'
    self.modelDescription = 'Fit SEIIR with Intervention on Daily Data with positive StudentT likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 8
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']), 
            ('Z', *self.defaults['Z']), 
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']),
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('cdof', *self.defaults['cdof'])
            )
    
    return js
