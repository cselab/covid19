#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   16/3/2020
# Email:  garampat@ethz.ch
# Description: Compute credible intervals

import numpy as np
from scipy.stats import norm
from scipy import optimize



def compute_credible_intervals( Z, Std, r):
  meanZ = np.mean( Z, axis=0 )
  Ns = Z.shape[0]
  n  = Z.shape[1]

  iv = {}
  iv['Intervals'] = []
  # for each time point in the `data` solve a non-linear equation to find the credible interval
  for p in r:

    local_iv = {}
    local_iv['Percentage'] = p
    local_iv['Low Interval']  = []
    local_iv['High Interval'] = []

    p_low  = 0.5-p/2
    p_high = 0.5+p/2

    for k in range(n):
      if( np.any( Std[:,k]==0. ) ):
        # local_iv['Low Interval'].append( np.nan )
        # local_iv['High Interval'].append( np.nan )
        local_iv['Low Interval'].append( meanZ[k] )
        local_iv['High Interval'].append( meanZ[k] )
      else:
        F = lambda x: np.mean( norm.cdf( x, Z[:,k], Std[:,k] ) ) - p_low
        local_iv['Low Interval'].append( optimize.fsolve( F, meanZ[k] )[0] )
        F = lambda x: np.mean( norm.cdf( x, Z[:,k], Std[:,k] ) ) - p_high
        local_iv['High Interval'].append( optimize.fsolve( F, meanZ[k] )[0] )

    iv['Intervals'].append( local_iv)

  iv['Median'] = []
  for k in range(n):
    if( np.any( Std[:,k]==0. ) ):
      # iv['Median'].append( np.nan )
      iv['Median'].append( meanZ[k] )
    else:
      F = lambda x: np.sum( norm.cdf( x, Z[:,k], Std[:,k] ) ) / Ns - 0.5
      iv['Median'].append( optimize.fsolve( F, meanZ[k] )[0] )

  iv['Mean'] = meanZ.tolist()

  return iv
