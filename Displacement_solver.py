from scipy.optimize import fsolve
import numpy as np
from Bearing import Bearing
from Load import BearingLoads
from scipy import integrate

class Stiffness():
    def __init__(self, bearing: Bearing, alpha:float):
        self.bearing = bearing
        self.alpha = alpha
        self._gamma = self.gamma(self.alpha)
        self._Rx_i = self.Rx(self.bearing.inner, self.alpha)
        self._Rx_e = self.Rx(self.bearing.outer, self.alpha)
        self._Ry_i = self.Ry(self.bearing.inner)
        self._Ry_e = self.Ry(self.bearing.outer)
        self._R_i = 1/self.inverse_R(self._Rx_i, self._Ry_i)
        self._R_e = 1/self.inverse_R(self._Rx_e, self._Ry_e)
        self._GAMMA_i = self.GAMMA(self._Rx_i, self._Ry_i)
        self._GAMMA_e = self.GAMMA(self._Rx_e, self._Ry_e)
        self._kappa_i = self.solve_kappa(self._GAMMA_i)
        self._kappa_e = self.solve_kappa(self._GAMMA_e)
        self.K_i = self.K(self._kappa_i, self._R_i)
        self.K_e = self.K(self._kappa_e, self._R_e)
        self.K_n = self.Kn()

    def gamma(self, alpha):
        return self.bearing.ball.D * np.cos(alpha)/self.bearing._dm
    
    def Rx(self, part, alpha):
        #Verified
        if part.type =='Inner':
            Rx = (1-self.gamma(alpha)) * self.bearing.ball.D/2
        elif part.type == 'Outer':
            Rx = (1+self.gamma(alpha)) * self.bearing.ball.D/2
        return Rx
    
    def Ry(self, part):
        #Verified
        if part.type =='Inner':
            Ry = self.bearing.fi*self.bearing.ball.D/(2*self.bearing.fi -1)
        elif part.type == 'Outer':
            Ry = self.bearing.fo*self.bearing.ball.D/(2*self.bearing.fo -1)
        return Ry
    
    def GAMMA(self, Rx, Ry):
        """Curvature difference computation
        Verified

        Args:
            Rx (float): X curvature radius (shaft axis)
            Ry (float): Y curvature radius (ball rolling direction)

        Returns:
            float: Major gamma function ie curvature difference
        """
        return (1/self.inverse_R(Rx, Ry)) * ( 1/Rx - 1/Ry)
    
    def inverse_R(self, Rx, Ry):
        #Verified
        return 1/Rx + 1/Ry
    
    def F(self, k):
        
        """F(k) function. Elliptic integral of first kind
        Verified

        Args:
            k (float): Ellipse elongation

        Returns:
            float: F(k)
        """
        integrand = lambda t: 1/(np.sqrt(1-(1-1/k**2)*np.sin(t)**2))
        res, abs = integrate.quad(integrand, 0, np.pi/2, limit=200)
        return res

    def S(self, k):
        """S(k) function. Elliptic integral of first kind
        Verified

        Args:
            k (float): Ellipse elongation

        Returns:
            float: S(k)
        """
        integrand = lambda t: np.sqrt(1 - (1 - (1/k)**2)*np.sin(t)**2)
        res, abs = integrate.quad(integrand, 0, np.pi/2, limit=200)
        return res
    
    def H1(self, k, GAMMA):
        """Function to solve for k 
        Verified
        Args:
            k (float): Ellipse elongation
            GAMMA (float): Curvature difference [1/m]

        Returns:
            float: Gap with zero
        """                
        return 1 - 2/(k**2-1) * ((self.F(k)/self.S(k))-1) - GAMMA
    
    def solve_kappa(self, GAMMA):
        """Résout H1(k, GAMMA) = 0 pour k avec fsolve.
        Verified
        """
        try:
            res, infodict, ier, mesg = fsolve(
                self.H1,
                x0=2,  # Valeur initiale
                args=(GAMMA,),
                full_output=True
            )
            if ier == 1:  # Succès
                # print(f"Solution pour k: {res[0]}")
                return res[0]
            else:
                print(f"fsolve n'a pas convergé: {mesg}")
                return None
        except Exception as e:
            print(f"Erreur lors de la résolution: {e}")
            return None
        
    def equiv_E(self, part1, part2):
        """Equivalent modulus of elasticity
        Verified
        Args:
            part1 (Ball or Raceway): First part
            part2 (Ball or Raceway): Second part

        Returns:
            float: Equivalent modulus of elasticity [MPa]
        """
        return 2/( (1-part1.Mat.Nu**2)/part1.Mat.E + (1-part2.Mat.Nu**2)/part2.Mat.E)
    
    def K(self, kappa, R):
        return np.pi/3 * kappa * self.equiv_E(self.bearing.ball, self.bearing.outer) * np.sqrt(2*self.S(kappa) * R / self.F(kappa)**3)

    def Kn(self):
        """Computation of Kn, the ball stiffness

        Returns:
            float: Ball stiffness [N / m^(3/2)]
        """
        # print('GAMMAi', GAMMAi)
        # print('GAMMAe', GAMMAe)
        # print('Exact Ki', Ki)
        # print('Exact Ke', Ke)
        Kn = (1/(self.K_i**(2/3)) + 1/(self.K_i**(2/3)))**(-3/2)
        return Kn

class disp_solve_static():
    def __init__(self, bearing: Bearing, loads : BearingLoads):
        self.__bearing = bearing
        self.__psi_array = self.__bearing.psi_array
        self.__Z = self.__bearing.Z
        self.__A = self.__bearing.A
        self.__alpha0 = self.__bearing.alpha0
        self.__Ri = self.__bearing._Ri
        self._Stiffness = Stiffness(self.__bearing, self.__alpha0)
        self.__Kn = self._Stiffness.K_n
        self.__dm = self.__bearing._dm
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
        raise Exception("\n No solution has been found")
        return None

class disp_solve_dynamic():
    def __init__(self, bearing: Bearing, loads : BearingLoads):
        self.__bearing = bearing
        self.__D = self.__bearing.ball.D
        self.__fi = self.__bearing.fi
        self.__fo = self.__bearing.fo
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
        self.__wi = loads.wi
        self.__wo = loads.wo
    
    def Mg(self):
        return (2*np.pi)**2/60**3 * self.__bearing.ball.Mat.Rho*np.pi*self.__bearing.ball.D**5 * nR * nm*np.sin(beta)

    def A1j(self, deltaA, theta, psi):
        return self.__bearing.B*self.__bearing.ball.D*np.sin(self.__bearing.alpha0) + deltaA + theta*self.__bearing.Ri*np.cos(psi)
    
    def A2j(self, deltaR, psi):
        return self.__bearing.B*self.__bearing.ball.D*np.cos(self.__bearing.alpha0) + deltaR*np.cos(psi)

    def cos_alphaoj(self, X2j, deltaoj):
        return X2j/((self.__bearing.fo - 0.5)*self.__bearing.ball.D + deltaoj)
    
    def sin_alphaoj(self, X1j, deltaoj):
        return X1j/((self.__bearing.fo - 0.5)*self.__bearing.ball.D + deltaoj)
    
    def cos_alphaij(self, A2j, X2j, deltaij):
        return (A2j-X2j)/((self.__bearing.fi - 0.5)*self.__bearing.ball.D + deltaij)
    
    def sin_alphaij(self, A1j, X1j, deltaij):
        return (A1j-X1j)/((self.__bearing.fi - 0.5)*self.__bearing.ball.D + deltaij)
    
    def wm_wj(self):
        return
    
    def wr_wj(self):
        return
    
    def Mgj(self):
        return self.ball.J * self.wr_wj() * self.loads.w**2 * self.wm_wj() * np.sin(beta)
    
    def Fcj(self):
        return 0.5 * self.__bearing.ball.mass() * self.__bearing.dm * self.loads.w**2 * self.wm_wj()**2
    
    def solve_X_delta(self, x, A1j, A2j, psi):
        lambdaoj = 2 #For outer raceway control
        lambdaij = 0 #For outer raceway control
        X1j, X2j, deltaij, deltaoj = x[0], x[1], x[2], x[3]
        eq1 = (A1j - X1j)**2 + (A2j - X2j)**2 - ((self.__fi - 0.5)*self.__D + deltaij)**2
        eq2 = X1j**2 + X2j**2 - ((self.__fo - 0.5)*self.__D + deltaoj)**2
        eq3 = ((lambdaoj*self.Mgj()*X2j/self.__D) - self.Koj*deltaoj**1.5*X1j)/((self.__fo-0.5)*self.__D+deltaoj) + (self.Kij*deltaij**1.5*(A1j-X1j)-lambdaij*self.Mgj()/self.__D*(A2j-X2j))/((self.__fi-0.5)*self.__D+deltaij)
        eq4 = ((lambdaoj*self.Mgj()*X2j/self.__D) - self.Koj*deltaoj**1.5*X1j)/((self.__fo-0.5)*self.__D+deltaoj) + (self.Kij*deltaij**1.5*(A2j-X2j)-lambdaij*self.Mgj()/self.__D*(A1j-X1j))/((self.__fi-0.5)*self.__D+deltaij) - self.Fcj()

        return [eq1, eq2, eq3, eq4]

# Ceramic = Material('Steel')
# AMS5898 = Material('Steel')
# ball = Ball(Ceramic, 12.7/1000)
# inner_raceway = Raceway(AMS5898, 'Inner', 6.604/1000, 52.291/1000, 20/1000, b=20/1000, ds=60/1000)
# outer_raceway = Raceway(AMS5898, 'Outer', 6.604/1000, 77.706/1000, 70/1000, b=20/1000, ds=70/1000)
# bearing = Bearing(outer_raceway, inner_raceway, ball, Z=9, alpha0=0)
# loads = BearingLoads(500, 1000, 10)
# disp_solve(bearing, loads).solve_extended()
