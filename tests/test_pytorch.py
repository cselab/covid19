import sys,time 
import torch
import numpy as np
from scipy.integrate import solve_ivp
from decimal import Decimal
import matplotlib.pyplot as plt

sys.path.append('../')
import epidemics.ode_solver as solver

# ---------------------------------------------------------------------------- #
def sir_rhs_1(t, y, N, p ):
    S, I = y
    c1 = p[0] * S * I / N
    c2 = p[1] * I
    dSdt = -c1
    dIdt =  c1 - c2
    return dSdt, dIdt

def sir_rhs_2(t, y, N, p ):
    S, I = y
    c1 = p[0] * p[1] * S * I / N
    c2 = p[1] * I
    dSdt = -c1
    dIdt =  c1 - c2
    return dSdt, dIdt

# ---------------------------------------------------------------------------- #

# Sample parameters
y0 = [8476004,1]
p = [2.930196462548338, 5.25592088396661]
N = 8476005

t = np.arange(0,72,1)

t_1 = time.time()
sol_np = solver.solve_ode(sir_rhs_2,T=t[-1],y0=y0,args=(N,p),t_eval = t,backend='numpy',max_step=np.inf)
print(time.time()-t_1)

t_2 = time.time()
sol_t = solver.solve_ode(sir_rhs_2,T=t[-1],y0=y0,args=(N,p),t_eval = t,backend='torch',iterations_per_day=10)
print(sol_t.y[0])
print(time.time()-t_2)
fig = plt.figure()
plt.plot(sol_np.y[0],'-k',label='numpy')
plt.plot(sol_np.y[0],'--k',label='numpy')
plt.plot(sol_t.y[0].detach().numpy(),'-r',label='torch')
plt.plot(sol_t.y[0].detach().numpy(),'--r',label='torch')

# plt.plot(sol_np.y[0]-sol_t.y[0].detach().numpy())
plt.savefig('sir.pdf')

# ---------------------------------------------------------------------------- #
# SIR with intervention
def sir_rhs(t, y, N, p ):
    S, I = y

    if( t<p[3] ):
      beta = p[0]
    else:
      beta = p[0]*p[2]

    c1 = beta * p[1] * S * I / N
    c2 = p[1] * I
    dSdt = -c1
    dIdt =  c1 - c2
    return dSdt, dIdt

p = [1.1707981197206878, 2.198957884041561, 0.9087846462914353, 27.257789328279856, 1.4257123500712776]
y0 = [318713.0, 1.0]
N = y0[0]+y0[1]
t = np.arange(0,53,1)

sol_np = solver.solve_ode(sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = t,backend='numpy',max_step=1)
sol_t = solver.solve_ode(sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = t,backend='torch')

fig = plt.figure()
plt.plot(sol_np.y[0])
plt.plot(sol_t.y[0].detach().numpy())
plt.savefig('fig.pdf')

