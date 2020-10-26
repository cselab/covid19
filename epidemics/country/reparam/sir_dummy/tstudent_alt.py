import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sir_dummy.tstudent_alt'
    self.modelDescription = 'Fit SIR with interventions on Daily Infected Data with Positive StudentT likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )

 
  def get_variables_and_distributions( self ):
 
    self.nParameters = 2
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('cdof', *self.defaults['cdof'])
            )
    
    return js
