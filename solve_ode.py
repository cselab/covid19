#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   1/4/2020
# Email:  garampat@ethz.ch


import numpy as np
import matplotlib.pyplot as plt

# import epidemics.sir.modelBase as sir
# x = sir.modelBase(rawData=[1,2,3],populationSize=1000)
# N = int(1e3)
# T = 50
# sol = x.solve_ode( [N-100,100], T, N, [0.5,0.1] )
# t = np.linspace(0,T,1000)
# S = sol.sol(t)[0]
# I = sol.sol(t)[1]
#
# fig = plt.figure(figsize=(12, 8))
# ax  = fig.subplots(1)
# ax.plot( t, S,   '-', lw=2, label='S' )
# ax.plot( t, I,   '-', lw=2, label='I' )
# ax.grid()
# ax.legend(loc='upper left')
# # ax[k].set_yscale('log')
#
# plt.show()



import epidemics.seiir.modelBase as sir
x = sir.modelBase(rawData=[1,2,3],populationSize=1000)
N = int(1e3)
T = 30
sol = x.solve_ode( [N-100,0,100,0], T, N, [ 1.5, 0.5, 0.2, 3, 5 ] )
t = np.linspace(0,T,1000)
S = sol.sol(t)[0]
E = sol.sol(t)[1]
Ir = sol.sol(t)[2]
Iu = sol.sol(t)[3]

fig = plt.figure(figsize=(12, 8))
ax  = fig.subplots(1)
ax.plot( t, S,  '-', lw=2, label='S'  )
ax.plot( t, E,  '-', lw=2, label='E'  )
ax.plot( t, Ir, '-', lw=2, label='Ir' )
ax.plot( t, Iu, '-', lw=2, label='Iu' )



ax.grid()
ax.legend(loc='upper left')
# ax[k].set_yscale('log')

plt.show()
