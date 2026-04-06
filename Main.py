from Material import Material
from Ball import Ball
from Raceway import Raceway
from Bearing import Bearing
import numpy as np
import matplotlib.pyplot as plt
from Load import BearingLoads
from Displacement_solver import disp_solve
from Display import Display

Ceramic = Material('Steel')
AMS5898 = Material('Steel')
ball = Ball(Ceramic, 12.7/1000)
inner_raceway = Raceway(AMS5898, 'Inner', 6.604/1000, 52.291/1000, d_fit=40/1000, b=20/1000, ds=60/1000)
outer_raceway = Raceway(AMS5898, 'Outer', 6.604/1000, 77.706/1000, d_fit=90/1000, b=20/1000, ds=70/1000)
bearing = Bearing(outer_raceway, inner_raceway, ball, Z=9, alpha0=0)
loads = BearingLoads(500, 0, 0)

print('A:', bearing.A*1000,'mm')
print('Pd:', bearing.Pd*1000_000,'µm')
print('Pe:',bearing.Pe*1000,'mm')
disp=disp_solve(bearing, loads).solve_extended()
Display.Display(bearing)
Display.Display_ball_load(bearing, disp)
Display.Display_ball_pressure(bearing, disp)
Display.Display_ball_angle(bearing, disp)
plt.show()