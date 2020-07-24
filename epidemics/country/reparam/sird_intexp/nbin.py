import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.sird_intexp.nbin'
    self.modelDescription = 'Fit SIRD with exp intervention on Daily Data with Negative Binomial likelihood'
    self.likelihoodModel  = 'Negative Binomial'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 6
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']), 
            ('D', *self.defaults['D']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('k', *self.defaults['kexp']),
            ('r', *self.defaults['r']),
            )
    
    return js
