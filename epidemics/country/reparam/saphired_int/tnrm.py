import numpy as np
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.saphired_int.tnrm'
    self.modelDescription = 'Fit SAPHIRE on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 11
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Zl', *self.defaults['Zl']),
            ('Y', *self.defaults['Y']), 
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            ('Sigma', *self.defaults['Sigma'])
            )
    
    return js


