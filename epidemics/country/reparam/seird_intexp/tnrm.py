import numpy as np
from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seird_intexp.tnrm'
    self.modelDescription = 'Fit SEIRD on Daily Infected Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )


  def get_variables_and_distributions( self ):
 
    self.nParameters = 7
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('D', *self.defaults['D']),
            ('Z', *self.defaults['Z']), 
            ('eps', *self.defaults['eps']), 
            ('tact', *self.defaults['tact']),
            ('k', *self.defaults['k']),
            ('Sigma', *self.defaults['Sigma'])
            )
    
    return js


