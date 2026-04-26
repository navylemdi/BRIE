from Material import Material
from Ball import Ball
from Raceway import Raceway
from Bearing import Bearing
import numpy as np
import matplotlib.pyplot as plt
from Load import BearingLoads
from Displacement_solver import disp_solve_static, Stiffness
from Display import Display

Ceramic = Material('Steel')
AMS5898 = Material('Steel')
ball = Ball(Ceramic, 7.938/1000)
inner_raceway = Raceway(AMS5898, 'Inner', 4.25/1000, 40.546/1000, d_fit=35/1000, b=20/1000, ds=45/1000)
outer_raceway = Raceway(AMS5898, 'Outer', 4.29/1000, 56.448/1000, d_fit=62/1000, b=20/1000, ds=50/1000)
bearing = Bearing(outer_raceway, inner_raceway, ball, Z=10)
loads = BearingLoads(0, 2000, 0)





print('Alpha0:', np.degrees(bearing.alpha0), '°')
print('A:', bearing.A*1000,'mm')
print('Pd:', bearing._Pd*1_000_000,'µm')
print('Pe:', bearing._Pe*1000,'mm')
print('Kn:', Stiffness(bearing, bearing.alpha0).K_n, 'N/m1.5', Stiffness(bearing, bearing.alpha0).K_n*1000**-1.5, 'N/mm1.5')
print('Ke:', Stiffness(bearing, bearing.alpha0).K_e)
print('Ki:', Stiffness(bearing, bearing.alpha0).K_i)

# disp=disp_solve_static(bearing, loads).solve_extended()
# print('Force maximale:', np.round(bearing.Q_max(disp[0], disp[1], disp[2]),2), 'N' )
# print('Déplacement max outer:', np.round(bearing.deltaE(bearing.Q_max(disp[0], disp[1], disp[2]))*1000000,2), 'µm' )
# print('Déplacement max inner:', np.round(bearing.deltaI(bearing.Q_max(disp[0], disp[1], disp[2]))*1000000,2), 'µm' )
Display.Display(bearing)
# Display.Display_ball_load(bearing, disp)
# Display.Display_ball_pressure(bearing, disp)
# Display.Display_ball_angle(bearing, disp)
# Display.Display_contact_deformation(bearing, disp)
plt.show()