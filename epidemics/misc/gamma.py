#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from scipy.special import gamma
import scipy.integrate as integrate

def finiteDifference(x, eps=10e-6):
    return (gamma(x+eps)-gamma(x))/(eps)

def analytical(x, top=10e3):
    res = integrate.quad(lambda t: np.exp(-t)*(t**(x-1.))*np.log(t), 0, top)
    return res[0]
    

if __name__ == '__main__':

    N = 1000
    upper = 5
    print("Plotting Gamma function and its derivative on (0,{0}]".format(upper))

    x = np.linspace(0,upper,N)
    
    g = gamma(x)
    dg0 = [analytical(xi) for xi in x]
    dg1 = finiteDifference(x,1e-3)
    dg2 = finiteDifference(x,1e-5)
    dg3 = finiteDifference(x,1e-7)

    print("negative E-M constant and FD approximations:")
    print(-0.5772156649)
    print(finiteDifference(1, 1e-3))
    print(finiteDifference(1, 1e-5))
    print(finiteDifference(1, 1e-7))

    plt.plot(x,g)
    plt.show()

    plt.plot(x,dg1,x,dg2,x,dg3)
    plt.show()
    
