#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   1/4/2020
# Email:  garampat@ethz.ch
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

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



from epidemics.seiir.altone_tnrm import Model
x = Model()
N = x.regionalData.populationSize
T = x.regionalData.time[-1]
p = [ 2, 0.9, 0.05, 3, 5 ]

sol = x.solve_ode( [N-1,0,1,0], T, N, p )
t = np.linspace(0,T,1000)
S = sol.sol(t)[0]
E = sol.sol(t)[1]
Ir = sol.sol(t)[2]
Iu = sol.sol(t)[3]

fig = plt.figure(figsize=(12, 8))
ax  = fig.subplots(1)

ax.plot( t[1:], -p[2]*(np.diff(S)+np.diff(E)), '-', lw=2, label='Daily Reported' )



ax.plot( x.data['Model']['x-data'], x.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')



ax.grid()
ax.legend(loc='upper left')
plt.show()
# ax[k].set_yscale('log')

# plt.savefig('ode.pdf')
