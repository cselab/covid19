import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seird_intsmooth.geo'
    self.modelDescription = 'Fit SEIRD with Intervention on Daily Infected Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Geometric'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 7
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Z', *self.defaults['Z']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            )
    
    return js
