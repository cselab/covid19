#!/usr/bin/env python

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import time

def Solve(beta=2, Z=0.5, D=0.5, N=1000, I=1, tmax=30):
    beta = np.array(beta).flatten()
    Z = np.array(Z).flatten()
    D = np.array(D).flatten()
    N = np.array(N).flatten()
    I = np.array(I).flatten()

    def rhs(t, y):
        S, E, I = y.reshape((3, len(y) // 3))
        mu = 0
        return np.array([
                mu*(N-S)  -beta * I * S / N,
                beta * I * S / N - E * (1 / Z + mu),
                E / Z - I * (1 / D + mu),
                ]).flatten()
    y0 = np.array([N, np.zeros_like(N), I])
    sol = solve_ivp(rhs, [0, tmax], y0.flatten(), vectorized=True, max_step=0.1)

    return sol.t, (sol.y).reshape((y0.shape[0], y0.shape[1], len(sol.y[0])))

# russia
#data_y = [13,14,17,20,20,28,34,45,59,64,93,114,147,199,253,306,367,438,498,658,840,1036,1264,1534,1836,2337,2777]
# sourth korea
#data_y = [11,15,15,15,16,19,24,24,25,27,28,28,28,28,28,29,30,31,46,82,156,156,556,602,893,1146,1595,2022,2337,3526,4156,4812,5328,5766,6284,6767,7134,7382,7515,7515,7869,7979,8086,8162,8236,8320,8413,8413,8652,8799,8897,8961,9037,9137,9241,9332,9478,9583,9661,9786,9887,9976]

# parameters: beta, Z, D
p0 = [3, 3, 3]
# population
N0 = 8e6
# initial infected
I0 = 1
tmax = 60
# synthetic data
data_y = N0 - Solve(*p0, N0, I0, tmax=tmax)[1][0].flatten()
data_t = np.linspace(0, tmax, len(data_y))

nsmp = 1000
eps0 = 2
eps1 = 2
eps2 = 2
rbeta = np.array([p0[0]-eps0, p0[0]+eps0])
rZ = np.array([p0[1]-eps1, p0[1]+eps1])
rD = np.array([p0[2]-eps2, p0[2]+eps2])
rN = np.array([N0, N0])

vbeta = np.random.uniform(*rbeta, nsmp)
vZ = np.random.uniform(*rZ, nsmp)
vD = np.random.uniform(*rD, nsmp)
vN = np.random.uniform(*rN, nsmp)
vI = np.ones(nsmp) * I0


def rescale(a, r):
    if r[1] - r[0] == 0:
        return 0
    return np.clip((a - r[0]) / (r[1] - r[0])*1.6, 0., 1.)

time0 = time.time()

t, yy = Solve(vbeta, vZ, vD, vN, vI, tmax=tmax)

data_i = np.argmin(abs(t[:,np.newaxis] - np.array(data_t)[np.newaxis,:]), axis=0)


time1 = time.time()
print("solve: {:.2f} s".format(time1 - time0))

#plt.yscale('log')

for i in range(min(yy.shape[1], 10000)):
    color = [
            rescale(vbeta[i], rbeta),
            rescale(vZ[i], rZ),
            rescale(vD[i], rD),
            ]
    y = N0 - yy[0, i, :]
    #dist = ((y[data_i] - data_y) ** 2).mean() ** 0.5
    dist = np.mean(abs(y[data_i] - data_y)) / np.mean(data_y)
    like = np.clip(np.maximum(1., 1 - dist*3)**10, 0, 0.2)
    #like = (N0 - y[-1]) / N0
    plt.plot(t, y, alpha=like, lw=2, color=color)
    plt.title("""SEI model with N={:.0e},
beta={:}±{:}, Z={:}±{:}, D={:}±{:}""".format(
        N0, p0[0], eps0, p0[1], eps1, p0[2], eps2))
    plt.xlabel("$t$ [days]")
    plt.ylabel("$N - S$")

time2 = time.time()

#plt.ylim(0,3000)
plt.scatter(data_t, data_y, zorder=100, s=2, color='black')

print("plot: {:.2f} s".format(time2 - time1))

plt.savefig("sei.png")
