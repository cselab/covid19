#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   1/4/2020
# Email:  garampat@ethz.ch
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt


from epidemics.sir_d.model_base import ModelBase
x = ModelBase()
N = 1000
T = 10
p = [ 2, 1 ]
y0 = [N-1,1,0,0,0,0]
sol = x.solve_ode( y0, T, N, p )
t = np.linspace(0,T,1000)

fig = plt.figure(figsize=(12, 8))
ax  = fig.subplots(3)


S = sol.sol(t)[0]
I = sol.sol(t)[1]
dS0 = sol.sol(t)[2]
dI0 = sol.sol(t)[3]
dS1 = sol.sol(t)[4]
dI1 = sol.sol(t)[5]

ax[0].plot( t, S, '-', lw=2, label='S' )
ax[0].plot( t, I, '-', lw=2, label='I' )

ax[1].plot( t, dS0, '-', lw=2, label='dS-R0' )
ax[1].plot( t, dI0, '-', lw=2, label='dI-R0' )

ax[2].plot( t, dS1, '-', lw=2, label='dS-gamma' )
ax[2].plot( t, dI1, '-', lw=2, label='dI-gamma' )




epsilon = 1e-5
p = [ 2+epsilon, 1 ]
sol = x.solve_ode( y0, T, N, p )
Se = sol.sol(t)[0]
Ie = sol.sol(t)[1]
ax[1].plot( t, (Se-S)/epsilon, '--', lw=2, label='dS-R0 (FD)' )
ax[1].plot( t, (Ie-I)/epsilon, '--', lw=2, label='dI-R0 (FD)' )

p = [ 2, 1+epsilon ]
sol = x.solve_ode( y0, T, N, p )
Se = sol.sol(t)[0]
Ie = sol.sol(t)[1]
ax[2].plot( t, (Se-S)/epsilon, '--', lw=2, label='dS-gamma (FD)' )
ax[2].plot( t, (Ie-I)/epsilon, '--', lw=2, label='dI-gamma (FD)' )




for k in ax:
  k.grid()
  k.legend(loc='upper left')

plt.show()
