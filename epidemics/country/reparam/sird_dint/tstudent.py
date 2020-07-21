import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):

  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sird_dint.tstudent'
    self.modelDescription = 'Fit SIR with interventions on Daily Data with Positive StudentT Likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )

 
  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('D', *self.defaults['D']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('dof', *self.defaults['dof'])
            )
    
    return js
