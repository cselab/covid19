import numpy as np

from .model_base import ModelBase

class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.saphired_ints.tstudent_alt'
    self.modelDescription = 'Fit SAPHIRE with Intervention on Daily Infected Data with Positive StudentT likelihood'
    self.likelihoodModel  = 'Positive StudentT'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 10
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Zl', *self.defaults['Zl']),
            ('Y', *self.defaults['Y']), 
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('cdof', *self.defaults['cdof'])
            )
    
    return js
