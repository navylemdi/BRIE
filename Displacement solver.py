from scipy.optimize import fsolve
import itertools
import numpy as np
import matplotlib.pyplot as plt
from Bearing import Bearing
from Material import Material
from Ball import Ball
from Raceway import Raceway
from Load import BearingLoads

class disp_solve():
    def __init__(self, bearing: Bearing, loads : BearingLoads):
        self.bearing = bearing
        self.psi_array = self.bearing.psi_array
        self.Z = self.bearing.Z
        self.A = self.bearing.A
        self.alpha0 = self.bearing.alpha0
        self.Ri = self.bearing.Ri
        self.Kn = self.bearing.Kn
        self.dm=self.bearing.dm
        self.Fa = loads.Fa
        self.Fr = loads.Fr
        self.M = loads.M

    def denom(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return np.sqrt((np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))**2)

    def Ddenom_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)

    def Ddenom_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return np.cos(psi)*(np.cos(self.alpha0) + DeltaRbar*np.cos(psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)
    
    def Ddenom_DThetabar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return self.Ri * np.cos(psi)*(np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)
    
    def u1(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return ((self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)-1)**1.5) * (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))
    
    def Du1_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi)) + (self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1)**1.5

    def Du1_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))

    def Du1_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi)) + self.Ri*np.cos(psi)*(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1)**1.5
    
    def u2(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return ((self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)-1)**1.5) * (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi)
    
    def Du2_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi)

    def Du2_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi) + (np.cos(psi)**2)*(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5
    
    def Du2_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi)
    
    def Df1_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.Du1_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
    
    def Df1_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.Du1_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
        
    def Df1_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.Du1_Dthetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum

    def Df2_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.Du2_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u2(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
    
    def Df2_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.Du2_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u2(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
    
    def Df2_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.Du2_Dthetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u2(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
    
    def Df3_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += np.cos(psi)*(self.Du1_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
    
    def Df3_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += np.cos(psi)*(self.Du1_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.Kn * self.A**1.5 * sum
    
    def Df3_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += np.cos(psi)*(self.Du1_Dthetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -(self.dm/2)*self.Kn * self.A**1.5 * sum
    
    def Jacobian(self, X):
        DeltaAbar, DeltaRbar, Thetabar = X[0], X[1], X[2]
        return np.array([[self.Df1_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar), self.Df1_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar), self.Df1_Dthetabar(DeltaAbar, DeltaRbar, Thetabar)],
                         [self.Df2_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar), self.Df2_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar), self.Df2_Dthetabar(DeltaAbar, DeltaRbar, Thetabar)],
                         [self.Df3_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar), self.Df3_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar), self.Df3_Dthetabar(DeltaAbar, DeltaRbar, Thetabar)]])
    
    def f1(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5 * (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi)) / self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)
        return self.Fa - self.Kn * self.A**1.5 * sum
    
    def f2(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5 * (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi) / self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)
        return self.Fr - self.Kn * self.A**1.5 * sum
    
    def f3(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.psi_array:
            sum += (self.denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5 * (np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))*np.cos(psi) / self.denom(DeltaAbar, DeltaRbar, Thetabar, psi)
        return self.M - (self.dm/2)*self.Kn * self.A**1.5 * sum
    
    def system(self, X):
            DeltaAbar, DeltaRbar, Thetabar = X[0], X[1], X[2]
            return [self.f1(DeltaAbar, DeltaRbar, Thetabar), self.f2(DeltaAbar, DeltaRbar, Thetabar), self.f3(DeltaAbar, DeltaRbar, Thetabar)]
    
    def solve(self, x0):
        """Compute the displacement of inner ring

        Args:
            Fa (float): Axial force applied on the inner ring [N]
            Fr (float): Radial force applied on the inner ring [N]
            M (float): Moment applied on the inner ring [Nm]
        """
        # x0=[1, 0, 1]#[(self.Fa/self.Kn)**(1/1.5)/self.A, (self.Fr/self.Kn)**(1/1.5)/self.A, 1]
        res, infodict, ier, mesg = fsolve(self.system, x0, fprime=self.Jacobian, full_output=True)
        print('x0', x0)
        print('f(res)', self.system(res))
        print('result', res)
        print('result mm', res[0] * self.A *1000, res[1] * self.A*1000, np.degrees(res[2] * self.A))
        print(mesg)
        print(infodict)
        return res[0], res[1], res[2], ier
    
    def solve_extended(self):
        lim=5
        values = np.linspace(-5, 5, lim*2+1)
        for i in values:
            for j in values:
                for k in values:
                    solution = self.solve([i, j, k])
                    if solution[3] == 1:
                        return solution[:3]
        print("No solution has been found")
        return None



Ceramic = Material('Steel')
AMS5898 = Material('Steel')
ball = Ball(Ceramic, 12.7/1000)
inner_raceway = Raceway(AMS5898, 'Inner', 6.604/1000, 52.291/1000, 20/1000, b=20/1000, ds=60/1000)
outer_raceway = Raceway(AMS5898, 'Outer', 6.604/1000, 77.706/1000, 70/1000, b=20/1000, ds=70/1000)
bearing = Bearing(outer_raceway, inner_raceway, ball, Z=9, alpha0=0)
loads = BearingLoads(500, 1000, 10)
var = np.linspace(0, 2, 5)
disp_solve(bearing, loads).solve_extended()
# psi_array = solver.psi_array

# DeltaAbar = 2e-0
# DeltaRbar = 2e-0
# Thetabar = 0*2e-0
# psi = 3*np.pi/4

# cl=['r', 'g', 'b', 'y', 'm', 'k', 'c', 'grey', 'brown']
# linstyle=['solid', 'dotted', 'dashed', 'dashdot',(0, (3, 10, 1, 10, 1, 10))]
# denom_array=np.zeros((len(var), len(var), len(var), len(psi_array)))
# print(denom_array)
# for i1, DeltaAbar in enumerate(var):
#     for i2, DeltaRbar in enumerate(var):
#         for i3, Thetabar in enumerate(var):
#             for i4, psi in enumerate(psi_array):
#                 denom_array[i1, i2, i3, i4] = solver.denom(DeltaAbar, DeltaRbar, Thetabar, psi)
# plt.figure()
# for i2, DeltaRbar in enumerate(var):
#         for i3, Thetabar in enumerate(var):
#             for i4, psi in enumerate(psi_array):
#                 plt.plot(var, denom_array[:, i2, i3, i4], color=cl[i4], ls=linstyle[i3])
# plt.grid()
# plt.show()
# print(solver.denom(DeltaAbar, DeltaRbar, Thetabar, psi))
# print(solver.Jacobian(DeltaAbar, DeltaRbar, Thetabar))
# func = [solver.u2(DeltaAbar, DeltaRbar, i, psi) for i in var]
# func_deriv = [solver.Du2_Dthetabar(DeltaAbar, DeltaRbar, i, psi) for i in var]
# dx=var[1]-var[0]
# func_numerical_deriv = np.gradient(func, dx)
# print(func)

# plt.figure()
# plt.plot(var, func, label='Function')
# plt.plot(var, func_deriv, label='Derivation')
# plt.plot(var, func_numerical_deriv, label='Numerical derivation')
# plt.legend()
# plt.show()

