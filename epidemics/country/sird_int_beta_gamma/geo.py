import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.sird_int_beta_gamma.geo'
    self.modelDescription = 'Fit SIRD with interventions on Daily Data with Geometric Likelihood'
    self.likelihoodModel  = 'Geometric'

    super().__init__( **kwargs )

 
  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('beta', *self.defaults['beta']), 
            ('gamma', *self.defaults['gamma']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('dtact', *self.defaults['dtact']),
            ('kbeta', *self.defaults['kbeta']),
            )
    
    return js
