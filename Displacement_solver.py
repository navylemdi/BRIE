from scipy.optimize import fsolve
import numpy as np
from Bearing import Bearing
from Load import BearingLoads

class disp_solve():
    def __init__(self, bearing: Bearing, loads : BearingLoads):
        self.__bearing = bearing
        self.__psi_array = self.__bearing.psi_array
        self.__Z = self.__bearing.Z
        self.__A = self.__bearing.A
        self.__alpha0 = self.__bearing.alpha0
        self.__Ri = self.__bearing.Ri
        self.__Kn = self.__bearing.K_n
        self.__dm = self.__bearing.dm
        self.__Fa = loads.Fa
        self.__Fr = loads.Fr
        self.__M = loads.M

    def __denom(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return np.sqrt((np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi))**2 + (np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))**2)

    def __Ddenom_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)

    def __Ddenom_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return np.cos(psi)*(np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)
    
    def __Ddenom_DThetabar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return self.__Ri * np.cos(psi)*(np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)
    
    def __u1(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return ((self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)-1)**1.5) * (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi))
    
    def __Du1_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.__Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi)) + (self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1)**1.5

    def __Du1_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.__Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi))

    def __Du1_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.__Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi)) + self.__Ri*np.cos(psi)*(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1)**1.5
    
    def __u2(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return ((self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)-1)**1.5) * (np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi)
    
    def __Du2_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.__Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi)

    def __Du2_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.__Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi) + (np.cos(psi)**2)*(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5
    
    def __Du2_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        "Verified"
        return 1.5 * self.__Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * np.sqrt(self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - 1) * (np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi)
    
    def __Df1_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__Du1_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
    
    def __Df1_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__Du1_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
        
    def __Df1_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__Du1_Dthetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum

    def __Df2_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__Du2_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u2(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
    
    def __Df2_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__Du2_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u2(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
    
    def __Df2_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__Du2_Dthetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u2(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
    
    def __Df3_DdeltaAbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += np.cos(psi)*(self.__Du1_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
    
    def __Df3_DdeltaRbar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += np.cos(psi)*(self.__Du1_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -self.__Kn * self.__A**1.5 * sum
    
    def __Df3_Dthetabar(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += np.cos(psi)*(self.__Du1_Dthetabar(DeltaAbar, DeltaRbar, Thetabar, psi) * self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) - self.__u1(DeltaAbar, DeltaRbar, Thetabar, psi)*self.__Ddenom_DThetabar(DeltaAbar, DeltaRbar, Thetabar, psi))/self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)**2
        return -(self.__dm/2)*self.__Kn * self.__A**1.5 * sum
    
    def __Jacobian(self, X):
        DeltaAbar, DeltaRbar, Thetabar = X[0], X[1], X[2]
        return np.array([[self.__Df1_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar), self.__Df1_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar), self.__Df1_Dthetabar(DeltaAbar, DeltaRbar, Thetabar)],
                         [self.__Df2_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar), self.__Df2_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar), self.__Df2_Dthetabar(DeltaAbar, DeltaRbar, Thetabar)],
                         [self.__Df3_DdeltaAbar(DeltaAbar, DeltaRbar, Thetabar), self.__Df3_DdeltaRbar(DeltaAbar, DeltaRbar, Thetabar), self.__Df3_Dthetabar(DeltaAbar, DeltaRbar, Thetabar)]])
    
    def __f1(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5 * (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi)) / self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)
        return self.__Fa - self.__Kn * self.__A**1.5 * sum
    
    def __f2(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5 * (np.cos(self.__alpha0) + DeltaRbar*np.cos(psi))*np.cos(psi) / self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)
        return self.__Fr - self.__Kn * self.__A**1.5 * sum
    
    def __f3(self, DeltaAbar, DeltaRbar, Thetabar):
        sum=0
        for psi in self.__psi_array:
            sum += (self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi) -1)**1.5 * (np.sin(self.__alpha0) + DeltaAbar + self.__Ri*Thetabar*np.cos(psi))*np.cos(psi) / self.__denom(DeltaAbar, DeltaRbar, Thetabar, psi)
        return self.__M - (self.__dm/2)*self.__Kn * self.__A**1.5 * sum
    
    def __system(self, X):
            DeltaAbar, DeltaRbar, Thetabar = X[0], X[1], X[2]
            return [self.__f1(DeltaAbar, DeltaRbar, Thetabar), self.__f2(DeltaAbar, DeltaRbar, Thetabar), self.__f3(DeltaAbar, DeltaRbar, Thetabar)]
    
    def solve(self, x0):
        """Compute the displacement of inner ring

        Args:
            Fa (float): Axial force applied on the inner ring [N]
            Fr (float): Radial force applied on the inner ring [N]
            M (float): Moment applied on the inner ring [Nm]
        """
        res, infodict, ier, mesg = fsolve(self.__system, x0, fprime=self.__Jacobian, full_output=True)
        # print('x0', x0)
        # print('f(res)', self.system(res))
        # print('result', res)
        # print('result mm', res[0] * self.A *1000, res[1] * self.A*1000, np.degrees(res[2] * self.A))
        # print(mesg)
        # print(infodict)
        return res[0], res[1], res[2], ier
    
    def solve_extended(self):
        lim=5
        values = np.linspace(-lim, lim, lim*2+1)
        compt = 0
        for i in values:
            for j in values:
                for k in values:
                    solution = self.solve([i, j, k])
                    compt += 1
                    print('Displacement calculation in progress: ', np.round(compt/(lim*2+1)**3*100,2), '%', end="\r")
                    if solution[3] == 1:
                        print('', end="\n")
                        return solution[:3]
        print("No solution has been found")
        return None



# Ceramic = Material('Steel')
# AMS5898 = Material('Steel')
# ball = Ball(Ceramic, 12.7/1000)
# inner_raceway = Raceway(AMS5898, 'Inner', 6.604/1000, 52.291/1000, 20/1000, b=20/1000, ds=60/1000)
# outer_raceway = Raceway(AMS5898, 'Outer', 6.604/1000, 77.706/1000, 70/1000, b=20/1000, ds=70/1000)
# bearing = Bearing(outer_raceway, inner_raceway, ball, Z=9, alpha0=0)
# loads = BearingLoads(500, 1000, 10)
# disp_solve(bearing, loads).solve_extended()
