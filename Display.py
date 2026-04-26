import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle, Arc
from Bearing import Bearing
from Load import BearingLoads
from Displacement_solver import Stiffness
import numpy as np

class Display():
    def __init__(self):
        pass
    
    def Display_ball_load(bearing: Bearing, Displacements):
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        Q=[]
        Stiffness_obj = Stiffness(bearing, bearing.alpha0)
        psi=list(bearing.psi_range())
        for angle in psi:
            Q.append(bearing.Q(Displacements[0], Displacements[1], Displacements[2], angle, Stiffness_obj.K_n))
        psi.append(psi[0])
        Q.append(Q[0])
        ax.plot(psi, Q)
        ax.set_xticks(np.linspace(0, 2*np.pi, bearing.Z, endpoint=False)[:bearing.Z])
        plt.title('Force on ball [N]')
    
    def Display_ball_angle(bearing: Bearing, Displacements):
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        A=[]
        psi=list(bearing.psi_range())
        for angle in psi:
            A.append(np.degrees(np.arcsin(bearing.alpha(Displacements[0], Displacements[1], Displacements[2], angle))))
        psi.append(psi[0])
        A.append(A[0])
        ax.plot(psi, A)
        ax.set_xticks(np.linspace(0, 2*np.pi, bearing.Z, endpoint=False)[:bearing.Z])
        plt.title('Contact angle on ball [°]')
    
    def Display_contact_deformation(bearing: Bearing, Displacements):
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        DeltaI=[]
        DeltaE=[]
        Stiffness_obj = Stiffness(bearing, bearing.alpha0)
        psi=list(bearing.psi_range())
        for angle in psi:
            DeltaI.append(bearing.delta(bearing.Q(Displacements[0], Displacements[1], Displacements[2], angle, Stiffness_obj.K_n), Stiffness_obj.K_i) * 1_000_000)
            DeltaE.append(bearing.delta(bearing.Q(Displacements[0], Displacements[1], Displacements[2], angle, Stiffness_obj.K_n), Stiffness_obj.K_e) * 1_000_000)
        psi.append(psi[0])
        DeltaI.append(DeltaI[0])
        DeltaE.append(DeltaE[0])
        ax.plot(psi, DeltaI, color='r', label='Inner race')
        ax.plot(psi, DeltaE, color='b', label='Outer race')
        ax.set_xticks(np.linspace(0, 2*np.pi, bearing.Z, endpoint=False)[:bearing.Z])
        plt.title('Contact deformation on ball [µm]')
        plt.legend()
    
    def Display_ball_pressure(bearing: Bearing, Displacements):
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        Stiffness_obj = Stiffness(bearing, bearing.alpha0)
        kappa_i = Stiffness_obj._kappa_i
        kappa_e = Stiffness_obj._kappa_e
        P_i=[]
        P_e=[]
        psi=list(bearing.psi_range())
        for angle in psi:
            Q = bearing.Q(Displacements[0], Displacements[1], Displacements[2], angle, Stiffness_obj.K_n)
            a_i = bearing.a(Q, kappa_i, bearing.inner)
            b_i = bearing.b(Q, kappa_i, bearing.inner)
            P_i.append(bearing.P(Q,a_i,b_i)*1e-6)
            a_e = bearing.a(Q, kappa_e, bearing.outer)
            b_e = bearing.b(Q, kappa_e, bearing.outer)
            P_e.append(bearing.P(Q,a_e,b_e)*1e-6)
        psi.append(psi[0])
        P_i.append(P_i[0])
        P_e.append(P_e[0])
        ax.plot(psi, P_i, color='r', label='Inner race')
        ax.plot(psi, P_e, color='b', label='Outer race')
        ax.set_xticks(np.linspace(0, 2*np.pi, bearing.Z, endpoint=False)[:bearing.Z])
        plt.title('Contact pressure on ball[MPa]')
        plt.legend()

    def Display(bearing: Bearing):
        def sagittas(part):
            return part.r + part.ds/2 - bearing._Ri
    
        def arc_cord(part):
            return 2*np.sqrt( part.r**2 - (-part.r + sagittas(part))**2 )
    
        def arc_opening_angle(part):
            return np.degrees( 2*np.arcsin(arc_cord(part)/(2*part.r)) )
        
        fig = plt.figure()
        ax1 = plt.subplot(2,1,1)
        # X normal view
        ## Ball Drawing
        for i in range(0, bearing.Z):
            c=Circle((bearing._dm/2 * np.cos(i*2*np.pi/bearing.Z), bearing._dm/2 * np.sin(i*2*np.pi/bearing.Z)), bearing.ball.D/2, color='k')
            ax1.add_patch(c)
        
        ## Inner ring drawing
        ax1.add_patch(Circle((0,0), bearing.inner.d_fit/2, ec='b', fc='w'))

        ## Outer ring drawing
        ax1.add_patch(Circle((0,0), bearing.outer.d_fit/2, ec='r', fill = False))

        ax1.axis('equal')
        ax1.legend(handles = [Line2D([0], [0], marker ='o', color='w', markerfacecolor='k', markersize=15, label='Ball'),
                              Line2D([0], [0], color='b', lw=2, label='Inner ring'),
                              Line2D([0], [0], color='r', lw=2, label='Outer ring')])
        
        ax2 = plt.subplot(2,1,2)
        lw=1

        # Cut view
        ## Ball Drawing
        c=Circle((0, bearing._dm/2), bearing.ball.D/2, color='k')
        ax2.add_patch(c)

        ## Inner ring drawing
        ax2.plot([-bearing.inner.b/2, bearing.inner.b/2], [bearing.inner.d_fit/2,bearing.inner.d_fit/2], linewidth=lw, color='b')
        ax2.plot([-bearing.inner.b/2, -bearing.inner.b/2], [bearing.inner.d_fit/2, bearing.inner.ds/2], linewidth=lw, color='b')
        ax2.plot([bearing.inner.b/2, bearing.inner.b/2], [bearing.inner.d_fit/2, bearing.inner.ds/2], linewidth=lw, color='b')
        ax2.plot([-bearing.inner.b/2, - arc_cord(bearing.inner)/2], [bearing.inner.ds/2, bearing.inner.ds/2], linewidth=lw, color='b')
        ax2.plot([arc_cord(bearing.inner)/2, bearing.inner.b/2], [bearing.inner.ds/2, bearing.inner.ds/2], linewidth=lw, color='b')
        ax2.add_patch(Arc((0, bearing._Ri), 2*bearing.inner.r, 2*bearing.inner.r, color='b', theta1=(3/2*180 - arc_opening_angle(bearing.inner)/2), theta2=1/2 * (3*180+arc_opening_angle(bearing.inner)), linewidth=lw))
        ax2.plot(0, bearing._Ri, marker='+', c='b')        
        
        ## Outer ring drawing
        ax2.plot([-bearing.outer.b/2, bearing.outer.b/2], [bearing.outer.d_fit/2,bearing.outer.d_fit/2], linewidth=lw, color='r')
        ax2.plot([-bearing.outer.b/2, -bearing.outer.b/2], [bearing.outer.d_fit/2, bearing.outer.ds/2], linewidth=lw, color='r')
        dsed=2* (bearing.inner.d/2 + bearing.ball.D - 60e-6) #Diamètre d'épaulement à droite (Hauteur du passage pour une bille avec 10µm d'interference)
        ax2.plot([bearing.outer.b/2, bearing.outer.b/2], [bearing.outer.d_fit/2, dsed/2], linewidth=lw, color='r')
        theta1 = np.degrees(np.arcsin((dsed-2*bearing._Ro)/(2*bearing.outer.r)))
        theta2 = np.degrees(np.pi - np.arcsin((bearing.outer.ds-2*bearing._Ro)/(2*bearing.outer.r)))
        ax2.add_patch(Arc((0, bearing._Ro), 2*bearing.outer.r, 2*bearing.outer.r, color='r', theta1=theta1, theta2=theta2, linewidth=lw))
        ax2.plot([-bearing.outer.b/2, -bearing.outer.r*np.cos(np.pi-np.radians(theta2))], [bearing.outer.ds/2, bearing.outer.ds/2], lw=lw, color='r')
        ax2.plot([bearing.outer.b/2, bearing.outer.r*np.cos(np.radians(theta1))], [dsed/2, dsed/2], lw=lw, color='r')

        ax2.plot(0, bearing._Ro, marker='+', c='r')

        ax2.legend(handles = [Line2D([0], [0], marker ='o', color='w', markerfacecolor='k', markersize=15, label='Ball'),
                              Line2D([0], [0], color='b', lw=2, label='Inner ring'),
                              Line2D([0], [0], color='r', lw=2, label='Outer ring')])
        ax2.axis('equal')