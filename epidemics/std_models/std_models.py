#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   26/3/2020
# Email:  garampat@ethz.ch
import numpy as np



def standardDeviationModelConst( sigma, t, d=[] ):
  return [ sigma for s in t];

def standardDeviationModelLinear( sigma, t, d=[] ):
  return [ sigma*s for s in t];

def standardDeviationModelSqrt( sigma, t, d=[] ):
  return [ np.sqrt(sigma)*s for s in t];

def standardDeviationModelProportionalData( sigma, t, d ):
  return [ np.sqrt(sigma)*s for s in d];

def standardDeviationModelError( sigma, t, d=[] ):
  sys.exit('Unknown Standard Deviation Model.')


standard_deviation_models = {
  0: standardDeviationModelConst,
  1: standardDeviationModelLinear,
  2: standardDeviationModelSqrt,
  3: standardDeviationModelProportionalData
}
