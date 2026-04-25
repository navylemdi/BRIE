import numpy as np
from Ball import Ball 
from Raceway import Raceway


class Bearing():
    """Bearing object class
    """
    def __init__(self, outer : Raceway, inner: Raceway, ball:Ball, Z:int):
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
        self._dm=self.dm()
        self.fi=self.conformity(self.inner)
        self.fo=self.conformity(self.outer)
        # print('RBx', (self.dm-self.ball.D*np.cos(alpha0))/(2*np.cos(self.alpha0)))
        self.B = self.fi+self.fo-1
        self.A = self.B * self.ball.D
        self._Pd = self.Pd()
        self.alpha0 = self.alpha0_calc()
        self._Pe = self.Pe()
        
        self.psi_array = self.psi_range()
        self._Ri = self.Ri()
        self._Ro = self.Ro()

    def conformity(self, ring: Raceway):
        """Bearing conformity computation

        Args:
            ring (Raceway): One of both bearing ring (inner or outer)

        Returns:
            float: Conformity
        """
        return ring.r/self.ball.D
    
    def dm(self):
        return np.mean([self.outer.d, self.inner.d])
    
    def alpha0_calc(self):
        "Verified"
        return np.arccos(1-self._Pd/(2*self.A))
    
    def Pd(self):
        return self.outer.d - self.inner.d - 2*self.ball.D
        # return 2*self.B*self.ball.D * (1-np.cos(self.alpha0))
    
    def Pe(self):
        "Verified"
        return 2*self.A*np.sin(self.alpha0)
    
    def Ri(self):
        """Generate Radius of locus of inner raceway groove curvature center

        Returns:
            float: Radius of locus of raceway groove curvature centers
        """
        return self._dm/2 + (self.inner.r - self.ball.D/2) * np.cos(self.alpha0)
    
    def Ro(self):
        """Generate Radius of locus of outer raceway groove curvature center

        Returns:
            float: Radius of locus of raceway groove curvature centers
        """
        return self._Ri - self.A*np.cos(self.alpha0)
    
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
    
    def delta(self, Q, Stiffness):
        return (Q/Stiffness)**(1/1.5)
    
    def k_approx(self):
        """Compute the approximation of k according to GUAY
        Verified

        Returns:
            float: ki and ke
        """
        Rxi = self.Rx_i
        Ryi = self.Ry_i
        Rxe = self.Rx_e
        Rye = self.Ry_e
        rhoi = Ryi/Rxi
        rhoe = Rye/Rxe
        # print('rhoi', rhoi)
        # print('rhoe', rhoe)
        # k_approx_i = (Ryi/Rxi)**(2/np.pi)
        # k_approx_e = (Rye/Rxe)**(2/np.pi)
        k_approx_i = 1.18*(rhoi)**(0.598) - 0.19
        k_approx_e = 1.18*(rhoe)**(0.598) - 0.19
        # print('ki approximated', k_approx_i)
        # print('ke approximated', k_approx_e)
        return k_approx_i, k_approx_e

    def F_approx(self):
        """Approximation of F(k)
        Verified
        """
        Rxi = self.Rx_i
        Ryi = self.Ry_i
        Rxe = self.Rx_e
        Rye = self.Ry_e
        rhoi = Ryi/Rxi
        rhoe = Rye/Rxe
        F_ki_approx = np.pi/2 + (np.pi/2-1) * np.log(rhoi)
        F_ke_approx = np.pi/2 + (np.pi/2-1) * np.log(rhoe)
        print('F(ki) approximated', F_ki_approx)
        print('F(ke) approximated', F_ke_approx)

    def S_approx(self):
        """Approximation of S(k)
        Verified
        """
        Rxi = self.Rx_i
        Ryi = self.Ry_i
        Rxe = self.Rx_e
        Rye = self.Ry_e
        S_ki_approx= 1 + (np.pi/2-1) / (Ryi/Rxi)
        S_ke_approx= 1 + (np.pi/2-1) / (Rye/Rxe)
        print('S(ki) approximated', S_ki_approx)
        print('S(ke) approximated', S_ke_approx)

    def psi_range(self):
        psi_range=np.linspace(0, 2*np.pi, self.Z, False)
        return np.insert(psi_range[:len(psi_range)//2+1], 0,psi_range[len(psi_range)//2+1:]-2*np.pi)
        
    def Q_max(self, deltaAbar, deltaRbar, Thetabar, Kn):
        return Kn * self.A**1.5 * (((np.sin(self.alpha0) + deltaAbar + self.Ri*Thetabar)**2 + (np.cos(self.alpha0) + deltaRbar)**2)**0.5-1)**1.5
    
    def a(self, Q, kappa, part):
        """Compute ellipse axis
        Verified
        Args:
            Q (_type_): _description_
            kappa (_type_): _description_
            part (_type_): _description_

        Returns:
            float: _description_
        """
        if part.type=='Inner':
            R = self.R_i
        elif part.type=='Outer':
            R = self.R_e
        
        Sk = self.S(kappa)
        E = self.equiv_E(self.ball, part)  
        # print('R', R)
        # print('self.S(kappa)', self.S(kappa))
        # print('kappa', kappa)
        # print('Q', Q)
        # print('E', E)   
        a = (6*kappa**2 * Sk * Q * R/(np.pi*E))**(1/3)
        return a
    
    def b(self, Q, kappa, part):
        """Compute ellipse axis
        Verified
        Args:
            Q (_type_): _description_
            kappa (_type_): _description_
            part (_type_): _description_

        Returns:
            float: _description_
        """
        if part.type=='Inner':
            R = self.R_i

        elif part.type=='Outer':
            R = self.R_e
        
        Sk = self.S(kappa)
        E = self.equiv_E(self.ball, part)  
        b = (6*Sk*Q*R/(np.pi*kappa*self.equiv_E(self.ball, part)))**(1/3)
        return b
    
    def Q(self, DeltaAbar, DeltaRbar, Thetabar, psi, Kn):
        return Kn * self.A**1.5 * (((np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))**2)**0.5 -1)**1.5
    
    def P(self, Q, a, b):
        #Verified
        return 3*Q/(2*np.pi*a*b)
    
    def alpha(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return (np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi)) / ((np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar * np.cos(psi))**2)**0.5
    
    def alpha2(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return (np.cos(self.alpha0) + DeltaRbar*np.cos(psi)) / ((np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar * np.cos(psi))**2)**0.5
    
    def alpha3(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return (np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi)) / (np.cos(self.alpha0) + DeltaRbar * np.cos(psi))