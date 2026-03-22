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
        # print('RBx', (self.dm-self.ball.D*np.cos(alpha0))/(2*np.cos(self.alpha0)))
        self.B = self.fi+self.fo-1
        self.A=A
        self.gamma = self.gamma()
        # print('gamma', self.gamma)
        self.Rx_i = self.Rx(self.inner)
        self.Rx_e = self.Rx(self.outer)
        self.Ry_i = self.Ry(self.inner)
        self.Ry_e = self.Ry(self.outer)
        # print(self.Rx_i)
        # print(self.Rx_e)
        # print(self.Ry_i)
        # print(self.Ry_e)
        self.R_i = 1/self.inverse_R(self.Rx_i, self.Ry_i)
        self.R_e = 1/self.inverse_R(self.Rx_e, self.Ry_e)
        self.GAMMA_i = self.GAMMA(self.Rx_i, self.Ry_i)
        self.GAMMA_e = self.GAMMA(self.Rx_e, self.Ry_e)

        self.Kn = self.Kn()
        self.Ri = self.Ri()
        self.Ro = self.Ro()

    def conformity(self, ring: Raceway):
        """Bearing conformity computation

        Args:
            ring (Raceway): One of both bearing ring (inner or outer)

        Returns:
            float: Conformity
        """
        return ring.r/self.ball.D
    
    def Ri(self):
        return self.dm/2 + (self.inner.r - self.ball.D/2) * np.cos(self.alpha0)
    
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
    
    def gamma(self):
        return self.ball.D * np.cos(self.alpha0)/self.dm
    
    def Rx(self, part):
        #Verified
        if part.type =='Inner':
            Rx = (1-self.gamma) * self.ball.D/2
        elif part.type == 'Outer':
            Rx = (1+self.gamma) * self.ball.D/2
        return Rx
    
    def Ry(self, part):
        #Verified
        if part.type =='Inner':
            Ry = self.fi*self.ball.D/(2*self.fi -1)
        elif part.type == 'Outer':
            Ry = self.fo*self.ball.D/(2*self.fo -1)
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
        return 1 - 2/(k**2-1) * ((self.F(k)/self.S(k))-1) - GAMMA #np.sqrt((2*self.F(k)-self.S(k)*(1+GAMMA))/(self.S(k)*(1-GAMMA)))-k#
    
    def solve_k(self, GAMMA):
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
        gamma = self.gamma
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
        gamma = self.gamma
        Rxi = self.Rx_i
        Ryi = self.Ry_i
        Rxe = self.Rx_e
        Rye = self.Ry_e
        S_ki_approx= 1 + (np.pi/2-1) / (Ryi/Rxi)
        S_ke_approx= 1 + (np.pi/2-1) / (Rye/Rxe)
        print('S(ki) approximated', S_ki_approx)
        print('S(ke) approximated', S_ke_approx)

    def Kn(self):
        """Computation of Kn, the ball stiffness

        Returns:
            float: Ball stiffness [N / m^(3/2)]
        """
        # print('GAMMAi', GAMMAi)
        # print('GAMMAe', GAMMAe)
        k_res_i = self.solve_k(self.GAMMA_i)
        k_res_e = self.solve_k(self.GAMMA_e)
        # print('Exact ki', k_res_i)
        # print('Exact ke', k_res_e)
        Ke = np.pi/3 * k_res_e * self.equiv_E(self.ball, self.outer) * np.sqrt(2*self.S(k_res_e) * self.R_e / self.F(k_res_e)**3)
        Ki = np.pi/3 * k_res_i * self.equiv_E(self.ball, self.inner) * np.sqrt(2*self.S(k_res_i) * self.R_i / self.F(k_res_e)**3)
        # print('Exact Ki', Ki)
        # print('Exact Ke', Ke)
        Kn = (1/(Ki**(2/3)) + 1/(Ke**(2/3)))**(-3/2)
        return Kn
    
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

    def psi_range(self):
        psi_range=np.linspace(0, 2*np.pi, self.Z, False)
        return np.insert(psi_range[:len(psi_range)//2+1], 0,psi_range[len(psi_range)//2+1:]-2*np.pi)
    
    def func(self, x, Fa, Fr, M, psi_range):
            deltaAbar, deltaRbar, thetabar = x[0], x[1], x[2]
            eq1=Fa
            eq2=Fr
            eq3=M     
            for psi in psi_range:
                # if psi >=0 and psi<=np.pi:
                    print('Angle', np.round(np.degrees(psi),2), '°')                
                    denom = np.sqrt((np.sin(self.alpha0) + deltaAbar + self.Ri * thetabar * np.cos(psi))**2 + (np.cos(self.alpha0) + deltaRbar*np.cos(psi))**2)
                    print('denom', denom)
                    # print(deltaAbar, self.Ri * thetabar * np.abs(np.cos(psi)), deltaRbar*np.abs(np.cos(psi)))
                    num1 = (denom-1)**1.5 * (np.sin(self.alpha0) + deltaAbar + self.Ri*thetabar*np.cos(psi))
                    num2 = (denom-1)**1.5 * (np.cos(self.alpha0) + deltaRbar*np.abs(np.cos(psi)))*np.cos(psi)
                    num3 = (denom-1)**1.5 * (np.sin(self.alpha0) + deltaAbar + self.Ri*thetabar*np.cos(psi))*np.cos(psi)
                    # print(- self.Kn*self.A**1.5 * num1/denom, - self.Kn*self.A**1.5 * num2/denom, - self.dm/2 * self.Kn*self.A**1.5 * num3/denom )
                    eq1 += - self.Kn*self.A**1.5 * num1/denom
                    eq2 += - self.Kn*self.A**1.5 * num2/denom
                    eq3 += - self.dm/2 * self.Kn*self.A**1.5 * num3/denom
            print([eq1, eq2, eq3])
            return [eq1, eq2, eq3]
    
    def solve_disp(self, Fa, Fr, M):
        """Compute the displacement of inner ring

        Args:
            Fa (float): Axial force applied on the inner ring [N]
            Fr (float): Radial force applied on the inner ring [N]
            M (float): Moment applied on the inner ring [Nm]
        """
        psi_range = self.psi_range()
        x0=[1e-6, 1e-6, 1e-6]
        res, infodict, ier, mesg = fsolve(self.func, x0, args=(Fa, Fr, M, psi_range,), full_output=True)
        print('f(res)', self.func(res, Fa, Fr, M, psi_range))
        print('result', res)
        print(mesg)
        return res[0], res[1], res[2]
    
    def Q_max(self, deltaAbar, deltaRbar, Thetabar):
        return self.Kn * self.A**1.5 * (((np.sin(self.alpha0) + deltaAbar + self.Ri*Thetabar)**2 + (np.cos(self.alpha0) + deltaRbar)**2)**0.5-1)**1.5
    
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
        E=self.equiv_E(self.ball, part)  
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
        E=self.equiv_E(self.ball, part)  
        b = (6*Sk*Q*R/(np.pi*kappa*self.equiv_E(self.ball, part)))**(1/3)
        return b
    
    def Q(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return self.Kn * self.A**1.5 * (((np.sin(self.alpha0) + DeltaAbar + self.Ri*Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar*np.cos(psi))**2)**0.5 -1)**1.5
    
    def Pmax(self, Q, a, b):
        #Verified
        return 3*Q/(2*np.pi*a*b)
    
    def alpha(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return (np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi)) / ((np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar * np.cos(psi))**2)**0.5
    
    def alpha2(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return (np.cos(self.alpha0) + DeltaRbar*np.cos(psi)) / ((np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi))**2 + (np.cos(self.alpha0) + DeltaRbar * np.cos(psi))**2)**0.5
    
    def alpha3(self, DeltaAbar, DeltaRbar, Thetabar, psi):
        return (np.sin(self.alpha0) + DeltaAbar + self.Ri * Thetabar*np.cos(psi)) / (np.cos(self.alpha0) + DeltaRbar * np.cos(psi))

    def Display_ball_load(self, Fa, Fr, M):
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        disp = self.solve_disp(Fa, Fr, M)
        Q=[]
        psi=list(self.psi_range())
        for angle in psi:
            Q.append(self.Q(disp[0], disp[1], disp[2], angle))
        psi.append(psi[0])
        Q.append(Q[0])
        ax.plot(psi, Q)

    def Display(self):
        def sagittas(part):
            return part.r + part.ds/2 - self.Ri
    
        def arc_cord(part):
            return 2*np.sqrt( part.r**2 - (-part.r + sagittas(part))**2 )
    
        def arc_opening_angle(part):
            return np.degrees( 2*np.arcsin(arc_cord(part)/(2*part.r)) )
        
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
        lw=0.00001
        # Ball Drawing
        c=Circle((0, self.dm/2), self.ball.D/2, color='k')
        ax2.add_patch(c)

        # Inner ring drawing
        ax2.add_patch(Rectangle((-self.inner.b/2, self.inner.d/2 ), self.inner.b, lw, color='b'))
        ax2.add_patch(Rectangle((-self.inner.b/2, self.inner.d/2 ), lw, self.inner.ds/2 - self.inner.d/2, color='b'))
        ax2.add_patch(Rectangle((self.inner.b/2, self.inner.d/2 ), lw, self.inner.ds/2 - self.inner.d/2, color='b'))
        ax2.add_patch(Rectangle((-self.inner.b/2, self.inner.ds/2 ), self.inner.b/2 - arc_cord(self.inner)/2, lw, color='b'))
        ax2.add_patch(Rectangle((arc_cord(self.inner)/2, self.inner.ds/2 ), self.inner.b/2 - arc_cord(self.inner)/2, lw, color='b'))
        ax2.add_patch(Arc((0, self.Ri), 2*self.inner.r, 2*self.inner.r, color='b', theta1=(3/2*180 - arc_opening_angle(self.inner)/2), theta2=1/2 * (3*180+arc_opening_angle(self.inner)), linewidth=2))
        ax2.plot(0, self.Ri, marker='+', c='b')        
        
        # Outer ring drawing
        ax2.add_patch(Rectangle((-self.outer.b/2, self.outer.d/2 ), self.outer.b, lw, color='r'))
        
        ax2.add_patch(Arc((0, self.Ro), self.ball.D, self.ball.D, color='r', theta1=90, theta2=180, linewidth=lw*1000))
        ax2.plot(0, self.Ro, marker='+', c='r')

        ax2.legend(handles = [Line2D([0], [0], marker ='o', color='w', markerfacecolor='k', markersize=15, label='Ball'),
                              Line2D([0], [0], color='b', lw=2, label='Inner ring'),
                              Line2D([0], [0], color='r', lw=2, label='Outer ring')])
        ax2.axis('equal')
