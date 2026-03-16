import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import time
from scipy import integrate

a=0.0001
b=0.0005
x=0.001
y=0.001
z=.005

A=1
B=(a**2+b**2-x**2-y**2-z**2)
C=(a**2*b**2-x**2*b**2-y**2*a**2-z**2*a**2-z**2*b**2)
D=-z**2*a**2*b**2
p=100
s=time.time()
print(np.roots([A,B,C,D]))
print('Analytical method:', (time.time()-s)*1_000_000,'ns')
print(max(np.roots([A,B,C,D])))
u=max(np.roots([A,B,C,D]))

def Hu(u, a,b,x,y,z):
    return x**2/(a**2+u) + y**2/(b**2+u) + z**2/u -1

def solve_u(a,b,x,y,z):
    if (x/a)**2 + (y/b)**2 < 1:
        return 0
    else:
        try:
            results = optimize.root_scalar(
            Hu,
            args=(a,b,x,y,z),
            bracket=[1e-6, 1e3],   # intervalle large
            method='brentq'
            ).root
            return results
        except ValueError:
            print(a,b,x,y,z)

s=time.time()
print(solve_u(a,b,x,y,z))
print('Scipy method:',(time.time()-s)*1_000_000,'ns')

def P(p, a, b, u, x, y, z):
    integrand= lambda t: (1 - x**2/(a**2 + t) - y**2/(b**2 + t) - z**2/t) / np.sqrt((a**2+t)*(b**2+t)*t)
    result, abserr = integrate.quad(integrand, u, np.inf, limit=200)
    return result * 3*p/(16*np.pi)
s=time.time()
P(p,a,b,u,x,y,z)
print('P calculation time:', (time.time()-s)*1_000,'µs')

def convert_time(time):
    second=time
    minutes=second//60
    hour=minutes//60
    print("Elapsed time:", hour, 'h', minutes,'m', second, 's')

convert_time(60)