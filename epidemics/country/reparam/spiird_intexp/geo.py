import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.spiird_intexp.geo'
    self.modelDescription = 'Fit SPIIRD with Intervention on Daily Infected Data with Geometric likelihood'
    self.likelihoodModel  = 'Geometric'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 7
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Y', *self.defaults['Y']), 
            ('alpha', *self.defaults['alpha']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('k', *self.defaults['k']),
            )
    
    return js
