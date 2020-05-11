#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   1/4/2020
# Email:  garampat@ethz.ch
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt



# from epidemics.sir_d.model_base import ModelBase
# Ny = 2
# Np = 2
# x = ModelBase()
# N = 1000
# T = 10
# p = [ 2, 1 ]
# y0 = [N-1,1,0,0,0,0]


from epidemics.seiir_d.model_base import ModelBase
Ny = 4
Np = 5
x = ModelBase()
N = 1000
T = 10
p = [ 2, 0.1, 0.4, 2, 3 ]
y0 = np.zeros(((Np+1)*Ny,))
y0[0] = N-10
y0[1] = 10


#===============================================================================
# Common for all models
#===============================================================================

sol = x.solve_ode( y0, T, N, p )

t = np.linspace(0,T,1000)

fig = plt.figure(figsize=(12, 8))
ax  = fig.subplots(Np+1)

Y = sol.sol(t)

for i in range(Ny):
  ax[0].plot( t, Y[i], '-', lw=2, label='dY'+str(i) )

for k in range(Np):
  for i in range(Ny):
    ax[k+1].plot( t, Y[(k+1)*Ny+i], '-', lw=2, label='dY'+str(i+1)+'P'+str(k+1) )


epsilon = 1e-5

for k in range(Np):
  pe = p.copy(); pe[k] = pe[k] + epsilon
  sol = x.solve_ode( y0, T, N, pe )
  Ye = sol.sol(t)

  for i in range(Ny):
    ax[k+1].plot( t, (Ye[i]-Y[i])/epsilon, '--', lw=2, label='dY'+str(i+1)+'P'+str(k+1)+' (FD)' )


for k in ax:
  k.grid()
  k.legend(loc='upper left')

plt.show()
