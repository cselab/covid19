import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir_ints.tstudent_alt'
    self.modelDescription = 'Fit SIR with step intervention on Daily Data with Positive StudentT Likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )

 
  def get_variables_and_distributions( self ):
 
    self.nParameters = 5
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('D', *self.defaults['D']), 
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('cdof', *self.defaults['cdof'])
            )
    
    return js
