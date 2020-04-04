#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   26/3/2020
# Email:  garampat@ethz.ch
import numpy as np


# XXX Fix me

def standardDeviationModelConst( p, t, d=[] ):
  return [ p[2] for s in t];

def standardDeviationModelLinear( p, t, d=[] ):
  return [ p[2]*s for s in t];

def standardDeviationModelSqrt( p, t, d=[] ):
  return [ np.sqrt(p[2])*s for s in t];

def standardDeviationModelProportionalData( p, t, d ):
  return [ np.sqrt(p[2])*s for s in d];

def standardDeviationModelError( p, t, d=[] ):
  sys.exit('Unknown Standard Deviation Model.')


standard_deviation_models = {
  0: standardDeviationModelConst,
  1: standardDeviationModelLinear,
  2: standardDeviationModelSqrt,
  3: standardDeviationModelProportionalData
}
