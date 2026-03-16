import numpy as np
from Ball import Ball 
from Raceway import Raceway
from scipy.optimize import fsolve, root_scalar
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle, Arc

class Bearing():
    """Bearing object class
    """
    def __init__(self, outer : Raceway, inner: Raceway, ball:Ball, Z:int, alpha0:float, dm: float, A:float):
        """Initialization of the Bearing object

        Args:
            outer (Raceway): Outer ring of the bearing
            inner (Raceway): Inner ring of the bearing
            ball (Ball): Ball of the bearing (all balls are the same)
            Z (int) : Number of balls
            alpha0 (float): Free angle of contact [°]
            dm (float): Ball pitch diameter [m]
            A (float): Distance between raceway groove curvature centers [m]
        """

        self.outer = outer
        self.inner = inner
        self.ball = ball
        self.Z = Z
        self.alpha0=np.radians(alpha0)
        self.dm=dm
        self.fi=self.conformity(self.inner)
        self.fo=self.conformity(self.outer)
        self.B = self.fi+self.fo-1
        self.A=A
        self.Ri=self.Ri()
        self.Ro=self.Ro()
        self.Kn = self.Kn()

    def conformity(self, ring: Raceway):
        """Bearing conformity computation

        Args:
            ring (Raceway): One of both bearing ring (inner or outer)

        Returns:
            float: Conformity
        """
        return ring.r/self.ball.D
    
    def Ri(self):
        return self.dm/2 + (self.inner.r - self.ball.D/2)*np.cos(self.alpha0)
    
    def Ro(self):
        return self.Ri - self.A*np.cos(self.alpha0)
    
    def s(self, psi, deltaA, deltaR, theta):
        return np.sqrt((self.A*np.sin(self.alpha0 + deltaA + self.Ri * theta * np.cos(psi)))**2 + (self.A * np.cos(self.alpha0) + deltaR * np.cos(psi))**2)
    
    def variable_bar(self, variable):
        return variable/self.A
    
    def variable_bar_2_variable(self, variable_bar):
        return variable_bar*self.A
    
    def list_variable_bar_2_variable(self, list_variable_bar):
        new_list=[]

        for var in list_variable_bar:
            new_list.append(self.variable_bar_2_variable(var))
        return new_list

    def deltaN(self, psi, deltaA, deltaR, theta):
        deltaAbar=self.variable_bar(deltaA)
        deltaRbar=self.variable_bar(deltaR)
        thetabar=self.variable_bar(theta)
        psi = np.radians(psi)
        return self.A * ( np.sqrt((np.sin(self.alpha0) + deltaAbar + self.Ri * thetabar * np.cos(psi))**2 + (np.cos(self.alpha0) + deltaRbar * np.cos(psi))**2) - 1)
    
    def GAMMA(self, Rx, Ry):
        """Curvature difference computation according to GUAY

        Args:
            Rx (float): X curvature radius (ball rolling direction)
            Ry (float): Y curvature radius (shaft axis)

        Returns:
            float: _description_
        """
        return (1/Rx-1/Ry)/(1/Rx+1/Ry)
    
    def inverse_R(self,Rx, Ry):
        return 1/Rx + 1/Ry
    
    def F(self, k):
        """F(k) function. Elliptic integral of first kind

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

        Args:
            k (float): Ellipse elongation

        Returns:
            float: S(k)
        """
        integrand = lambda t: np.sqrt(1 - (1 - (1/k)**2)*np.sin(t)**2)
        res, abs = integrate.quad(integrand, 0, np.pi/2, limit=200)
        return res
    
    def H1(self, k,GAMMA):
        """Function to solve for k 

        Args:
            k (float): Ellipse elongation
            GAMMA (float): Curvature difference [1/m]

        Returns:
            float: Gap with zero
        """
        return 1 - 2/(k**2-1) * ((self.F(k)/self.S(k))-1) - GAMMA
    
    def solve_k(self, GAMMA):            
            # res=fsolve(func, [1])
            res = root_scalar(
                    self.H1,
                    args=(GAMMA),
                    bracket=[1e-3, 1e3],   # intervalle large
                    method='brentq'
                ).root
            return res
    
    def k_approx(self):
        """Compute the approximation of k according to GUAY

        Returns:
            float: ki and ke
        """
        gamma = self.ball.D * np.cos(self.alpha0)/self.dm
        Rxi = (1-gamma)*self.ball.D/2
        Ryi = self.fi * self.ball.D/(2*self.fi-1)
        Rxe = (1+gamma)*self.ball.D/2
        Rye = self.fo * self.ball.D/(2*self.fo-1)
        # k_approx_i = (Ryi/Rxi)**(2/np.pi)
        # k_approx_e = (Rye/Rxe)**(2/np.pi)
        k_approx_i = 1.18*(Ryi/Rxi)**(0.598) - 0.19
        k_approx_e = 1.18*(Rye/Rxe)**(0.598) - 0.19
        print('ki approximated', k_approx_i)
        print('ke approximated', k_approx_e)
        return k_approx_i, k_approx_e

    def F_approx(self):
        """Approximation of F(k)
        """
        gamma = self.ball.D * np.cos(self.alpha0)/self.dm
        Rxi = (1-gamma)*self.ball.D/2
        Ryi = self.fi * self.ball.D/(2*self.fi-1)
        Rxe = (1+gamma)*self.ball.D/2
        Rye = self.fo * self.ball.D/(2*self.fo-1)
        F_ki_approx=np.pi/2 + (np.pi/2-1) * np.log(Ryi/Rxi)
        F_ke_approx=np.pi/2 + (np.pi/2-1) * np.log(Rye/Rxe)
        print('F(ki) approximated', F_ki_approx)
        print('F(ke) approximated', F_ke_approx)

    def S_approx(self):
        """Approximation of S(k)
        """
        gamma = self.ball.D * np.cos(self.alpha0)/self.dm
        Rxi = (1-gamma)*self.ball.D/2
        Ryi = self.fi * self.ball.D/(2*self.fi-1)
        Rxe = (1+gamma)*self.ball.D/2
        Rye = self.fo * self.ball.D/(2*self.fo-1)
        S_ki_approx=1 + (np.pi/2-1) / (Ryi/Rxi)
        S_ke_approx=1 + (np.pi/2-1) / (Rye/Rxe)
        print('S(ki) approximated', S_ki_approx)
        print('S(ke) approximated', S_ke_approx)

    def Kn(self):
        """Computation of Kn, the ball stiffness

        Returns:
            float: Ball stiffness
        """
        gamma = self.ball.D * np.cos(self.alpha0)/self.dm
        Rxi = (1-gamma)*self.ball.D/2
        Ryi = self.fi * self.ball.D/(2*self.fi-1)
        Rxe = (1+gamma)*self.ball.D/2
        Rye = self.fo * self.ball.D/(2*self.fo-1)
        inverse_Ri = self.inverse_R(Rxi,Ryi)
        inverse_Re = 1/Rxe+1/Rye
        GAMMAi = self.GAMMA(Rxi, Ryi)
        GAMMAe = self.GAMMA(Rxe, Rye)
        # print('GAMMAi', GAMMAi)
        # print('GAMMAe', GAMMAe)
        k_res_i = self.solve_k(GAMMAi)
        k_res_e = self.solve_k(GAMMAe)
        # print('Exact ki', k_res_i)
        # print('Exact ke', k_res_e)
        Ke = np.pi/3 * k_res_e * self.equiv_E(self.ball, self.outer) * np.sqrt(2*self.S(k_res_e) * 1/inverse_Re / self.F(k_res_e)**3)
        Ki = np.pi/3 * k_res_i * self.equiv_E(self.ball, self.inner) * np.sqrt(2*self.S(k_res_i) * 1/inverse_Ri / self.F(k_res_e)**3)
        # print('Exact Ki', Ki)
        # print('Exact Ke', Ke)
        Kn = (1/(Ki**(2/3)) + 1/(Ke**(2/3)))**(-3/2)
        return Kn
    
    def equiv_E(self, part1, part2):
        """Equivalent modulus of elasticity

        Args:
            part1 (Ball or Raceway): First part
            part2 (Ball or Raceway): Second part

        Returns:
            float: Equivalent modulus of elasticity [MPa]
        """
        return 2/((1-part1.Mat.Nu**2)/part1.Mat.E + (1-part2.Mat.Nu**2)/part2.Mat.E)

    def solve_disp(self, Fa, Fr, M):
        """Compute the displacement of inner ring

        Args:
            Fa (float): Axial force applied on the inner ring [N]
            Fr (float): Radial force applied on the inner ring [N]
            M (float): Moment applied on the inner ring [Nm]
        """
        psi_range = np.linspace(-np.pi, np.pi, self.Z, False)
        # print('psi range', psi_range,'rad')
        # print(psi_range*180/np.pi, '°')

        def func(x):
            deltaAbar, deltaRbar, thetabar = x[0], x[1], x[2]
            eq1=Fa
            eq2=Fr
            eq3=M            
            for psi in psi_range:
                # print(psi)
                denom = (np.sqrt((np.sin(self.alpha0) + deltaAbar + self.Ri*thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + deltaRbar*np.cos(psi))**2))
                
                num1 = (np.sqrt(((np.sin(self.alpha0) + deltaAbar + self.Ri * thetabar * np.cos(psi))**2 + (np.cos(self.alpha0) + deltaRbar * np.cos(psi))**2))-1)**1.5 * (np.sin(self.alpha0) + deltaAbar + self.Ri*thetabar*np.cos(psi))
                num2 = (np.sqrt(((np.sin(self.alpha0) + deltaAbar + self.Ri * thetabar * np.cos(psi))**2 + (np.cos(self.alpha0) + deltaRbar * np.cos(psi))**2))-1)**1.5 * (np.cos(self.alpha0) + deltaRbar*np.cos(psi))*np.cos(psi)
                num3 = (np.sqrt(((np.sin(self.alpha0) + deltaAbar + self.Ri * thetabar * np.cos(psi))**2 + (np.cos(self.alpha0) + deltaRbar * np.cos(psi))**2))-1)**1.5 * (np.sin(self.alpha0) + deltaAbar + self.Ri*thetabar*np.cos(psi))*np.cos(psi)
                
                eq1 += - self.Kn*self.A**1.5 * num1/denom
                eq2 += - self.Kn*self.A**1.5 * num2/denom
                eq3 += - self.dm/2 * self.Kn*self.A**1.5 * num3/denom 
            
            return [eq1, eq2, eq3]
        
        x0=[1, 1 , 1]
        res, infodict, ier, mesg = fsolve(func, x0, full_output=True)
        print('f(xo)', func(x0))
        print('f(res)', func(res))
        print('result', res)
        print(mesg)
        return res[0], res[1], res[2]
    
    def Q_max(self, deltaAbar, deltaRbar, Thetabar):
        return self.Kn * self.A**1.5 * (((np.sin(self.alpha0) + deltaAbar + self.Ri*Thetabar)**2 + (np.cos(self.alpha0) + deltaRbar)**2)**0.5-1)**1.5

    def Display(self):
        fig = plt.figure()
        ax1 = plt.subplot(2,1,1)
        # Ball Drawing
        for i in range(0, self.Z):
            c=Circle((self.dm/2 * np.cos(i*2*np.pi/self.Z), self.dm/2 * np.sin(i*2*np.pi/self.Z)), self.ball.D/2, color='k')
            ax1.add_patch(c)
        
        # Inner ring drawing
        ax1.add_patch(Circle((0,0), self.inner.d/2, ec='b', fc='w'))

        # Outer ring drawing
        ax1.add_patch(Circle((0,0), self.outer.d/2, ec='r', fill = False))
        ax1.axis('equal')
        ax1.legend(handles = [Line2D([0], [0], marker ='o', color='w', markerfacecolor='k', markersize=15, label='Ball'),
                              Line2D([0], [0], color='b', lw=2, label='Inner ring'),
                              Line2D([0], [0], color='r', lw=2, label='Outer ring')])
        
        ax2 = plt.subplot(2,1,2)
        
        # Ball Drawing
        c=Circle((0, self.dm/2), self.ball.D/2, color='k')
        ax2.add_patch(c)

        # Inner ring drawing
        ax2.add_patch(Rectangle((-self.inner.b/2, self.inner.d/2 ), self.inner.b, self.inner.b/50, color='b'))
        ax2.add_patch(Arc((0, self.Ri), self.inner.r+10/1000, self.inner.r+10/1000, color='b', theta1=-90, theta2=0))

        # Outer ring drawing
        ax2.add_patch(Rectangle((-self.outer.b/2, self.outer.d/2 ), self.outer.b, self.outer.b/50, color='r'))
        ax2.add_patch(Arc((0, self.Ro), self.outer.r, self.outer.r, color='r', theta1=90, theta2=180))


        ax2.legend(handles = [Line2D([0], [0], marker ='o', color='w', markerfacecolor='k', markersize=15, label='Ball'),
                              Line2D([0], [0], color='b', lw=2, label='Inner ring'),
                              Line2D([0], [0], color='r', lw=2, label='Outer ring')])
        ax2.axis('equal')
        plt.show()
