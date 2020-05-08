from scipy.integrate import solve_ivp
import sys
import numpy as np
from decimal import Decimal
import torch

'''
# ------------------------------ ODE Solver ---------------------------------- #
 
    This script handles the ODE part of the model. Provided a RHS the ODE is 
    solved either using pytorch or numpy/scipy.

# ---------------------------------------------------------------------------- #
'''

def solve_ode(f,T,y0,args,t_eval,backend='numpy',iterations_per_day=10,max_step=np.inf):
    if backend == 'numpy':
        sol = solve_ivp(f,t_span=[0, T],y0=y0, args=args,t_eval=t_eval,max_step=max_step)
    elif backend == 'torch':
        sol = solve_ivp_torch(f,T,y0,args,t_eval=t_eval,iterations_per_day=iterations_per_day)
    else:
        sys.exit(f"\n[ODE SOLVER] Unknown backend, choose (numpy or torch)\n")
    return sol
  
def solve_ivp_torch(f,T,y0,args, t_eval,iterations_per_day=10):
    '''
        ODE solver using pytorch. ~100x slower than scipy 
        Using fixed time step RK4 integration
    '''
    T = int(T)
    t_eval = [int(i) for i in t_eval]
    # Here somehow fix dt automatically in a smart way using values of p
    dt = 1/iterations_per_day
    assert((Decimal(str(T)) % Decimal(str(dt)))==0.0) # weird fix?

    y = torch.autograd.Variable(torch.DoubleTensor(y0),requires_grad=True)
    p = torch.autograd.Variable(torch.DoubleTensor(args[1]),requires_grad=True)
    N = torch.autograd.Variable(torch.DoubleTensor([args[0]]),requires_grad=True)
    out = torch.autograd.Variable(torch.empty(len(y0),T+1,dtype=torch.double),requires_grad=False)
    
    out[:,0] = y.clone()
    t_out = np.arange(0,T+1,1)
    t = 0
    for day in range(T):
        for it in range(iterations_per_day):
            y, t = RK4_step(f,y,t,dt,N,p)
        out[:,day+1] = y.clone()

    # Return at t_eval
    idx = np.where(t_out  == t_eval)
    sol = Solution(out[:,idx],t_eval,p)

    return sol

def FE_step(f, y, t, dt, N, p):
    return y + dt * torch.cat(f(t,y,N,p)), t + dt

def RK4_step(f, y, t, dt, N, p):
    k1 = dt * torch.cat(f(t,y,N,p))
    k2 = dt * torch.cat(f(t + 0.5 * dt, y + 0.5 * k1, N, p))
    k3 = dt * torch.cat(f(t + 0.5 * dt, y + 0.5 * k2, N, p))
    k4 = dt * torch.cat(f(t+dt, y + k3,N,p))
    return y + k1/6 + k2/3 + k3/3 + k4/6, t + dt

# -------------------------------- Utils ------------------------------------- #

class Solution():
    def __init__(self,y,t,p):
        self.y = y[0]
        self.t = t
        self.p = p

def get_gradients(y,p):
    # Computes the Jacobian of the output w.r.t the model params
    # Returns a [n,m] matrix n:output dim m: param dim
    n = len(y)
    m = len(p)

    J = np.zeros([n,m])
    for i in range(n):
        v = np.ones(n)
        v[i] = 1
        v = torch.from_numpy(v).float()
        y.backward(v, retain_graph = True)
        J[i,:] = p.grad.numpy()

    return J

def to_list(y):
    if isinstance(y,list): 
        return y
    elif isinstance(y,torch.Tensor): 
        return y.detach().numpy().tolist()
    else:
        return y.tolist()

def check_zeros(y,tol):
    return [tol if x<tol else tol for x in y]

def append_zero(y):
    if isinstance(y,list): 
        return [0] + y
    elif isinstance(y,torch.Tensor): 
        return torch.cat([torch.zeros(1),y])
    else:
        return [0] + y.tolist()




