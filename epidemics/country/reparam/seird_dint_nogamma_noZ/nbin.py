import numpy as np
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seird_dint_nogamma_noZ.nbin'
    self.modelDescription = 'Fit SEIRD on Daily Data with Negative Binomial Likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 5
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('eps', *self.defaults['eps']),
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('r', *self.defaults['r'])
            )
    
    return js


