import numpy as np

from .model_base import ModelBase


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'country.reparam.seiird_dint_nogamma_noZ.tnrm'
    self.modelDescription = 'Fit SEIIRD with Interventions on Daily Data with Positive Normal Likelihood'
    self.likelihoodModel  = 'Positive Normal'

    super().__init__( **kwargs )

  def get_variables_and_distributions( self ):
 
    self.nParameters = 7
    js = self.get_uniform_priors(
            ('R0', *self.defaults['R0']),
            ('mu', *self.defaults['mu']), 
            ('alpha', *self.defaults['alpha']),
            ('eps', *self.defaults['eps']),
            ('tact', *self.defaults['tact']),
            ('kbeta', *self.defaults['kbeta']),
            ('Sigma', *self.defaults['Sigma'])
            )
    
    return js
