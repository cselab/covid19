import sys,time 
import torch
import numpy as np
from scipy.integrate import solve_ivp
from decimal import Decimal

sys.path.append('../')
import epidemics.ode_solver as solver

# ---------------------------------------------------------------------------- #
def sir_rhs(t, y, N, p ):
  S, I = y
  c1 = p[0] * S * I / N
  c2 = p[1] * I
  dSdt = -c1
  dIdt =  c1 - c2
  return dSdt, dIdt
# ---------------------------------------------------------------------------- #

# Sample parameters
y0 = [800000,23]
p = [10,10]
N = 800000+23

t = np.arange(0,50,1)

t_1 = time.time()
sol_np = solver.solve_ode(sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = t,backend='numpy')
print(time.time()-t_1)

t_2 = time.time()
sol_torch = solver.solve_ode(sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = t,backend='torch')
print(time.time()-t_2)

