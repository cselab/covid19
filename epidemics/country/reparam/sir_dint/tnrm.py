import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):

  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir_dint.tnrm'
    self.modelDescription = 'Fit SIR with interventions on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

 
  def get_variables_and_distributions( self ):
 
    self.nParameters = 5
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('D', *self.defaults['D']), 
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('Sigma', *self.defaults['Sigma'])
            )
    
    return js
