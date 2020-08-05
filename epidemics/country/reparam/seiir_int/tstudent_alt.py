import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiir_int.tstudent_alt'
    self.modelDescription = 'Fit SEIIR with Interventions on Daily Infected Data with Positive StudentT Likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )

  def get_variables_and_distributions( self ):
 
    self.nParameters = 9
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']), 
            ('Z', *self.defaults['Z']), 
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']),
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            ('cdof', *self.defaults['cdof']),
            )
    
    return js
