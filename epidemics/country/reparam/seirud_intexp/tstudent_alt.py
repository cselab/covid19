import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seirud_intexp.tstudent_alt'
    self.modelDescription = 'Fit SEIRUD with Intervention on Daily Infected Data with Positive StudentT likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 9
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Z', *self.defaults['Z']),
            ('Y', *self.defaults['Y']), 
            ('alpha', *self.defaults['alpha']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('k', *self.defaults['k']),
            ('cdof', *self.defaults['cdof'])
            )
    
    return js
