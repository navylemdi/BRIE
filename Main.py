from Material import Material
from Ball import Ball
from Raceway import Raceway
from Bearing import Bearing
import numpy as np


Ceramic = Material('Ceramic')
AMS5898 = Material('AMS5898')
ball = Ball(Ceramic, 6.35/1000)
inner_raceway = Raceway(AMS5898, 3.397/1000, 20/1000, 12/1000)
outer_raceway = Raceway(AMS5898, 3.429/1000, 42/1000, 12/1000)
bearing = Bearing(outer_raceway, inner_raceway, ball, 14, 15, 30.994/1000, 0.481854321427134/1000)
# disp = bearing.solve_disp(1000,1000,2)
# print(bearing.Q_max(disp[0], disp[1], disp[2]), 'N')
print('Ri:',bearing.Ri)
print('Ro:',bearing.Ro)
print('Rm:',bearing.dm/2)
print(bearing.Display())





