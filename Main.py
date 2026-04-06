from Material import Material
from Ball import Ball
from Raceway import Raceway
from Bearing import Bearing
import numpy as np
import matplotlib.pyplot as plt

Ceramic = Material('Steel')
AMS5898 = Material('Steel')
ball = Ball(Ceramic, 12.7/1000)
inner_raceway = Raceway(AMS5898, 'Inner', 6.604/1000, 52.291/1000, 20/1000, b=20/1000, ds=60/1000)
outer_raceway = Raceway(AMS5898, 'Outer', 6.604/1000, 77.706/1000, 70/1000, b=20/1000, ds=70/1000)
bearing = Bearing(outer_raceway, inner_raceway, ball, Z=9, alpha0=0)
print('A:', bearing.A*1000,'mm')
print('Pd:', bearing.Pd*1000_000,'µm')
print('Pe:',bearing.Pe*1000,'mm')
Fa=0
Fr=0
M=10
disp = bearing.solve_disp(Fa, Fr, M)
# disp_mm = bearing.variable_bar_2_variable(disp)
print('disp [mm]', bearing.list_variable_bar_2_variable(disp))
# Qmax = bearing.Q_max(disp[0], disp[1], disp[2])
# bearing.Display()
bearing.Display_ball_load(Fa, Fr, M)
bearing.Display_ball_pressure(Fa, Fr, M)
plt.show()
# print(Qmax, 'N')
# alpha0=[]
# alpha=[]
# alpha2=[]
# alpha3=[]
# psi=np.linspace(-np.pi, np.pi, bearing.Z , False)

# for angle in psi:
#     alpha0.append(bearing.alpha0)
#     alpha.append(np.arcsin((bearing.alpha(disp[0], disp[1], disp[2], angle))))
#     alpha2.append(np.arccos(bearing.alpha2(disp[0], disp[1], disp[2], angle)))
#     alpha3.append(np.arctan(bearing.alpha3(disp[0], disp[1], disp[2], angle)))
# print(np.degrees(alpha))
# print(np.degrees(alpha2))
# print(np.degrees(alpha3))
# print(np.degrees(np.arcsin(bearing.alpha(0,0,0,0))))
# sum=[]
# for i in range(0, len(psi)):
#     sum.append(alpha[i]**2 + alpha2[i]**2)
# print(sum)

# fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
# ax.plot(psi, np.degrees(alpha0), c='k')
# ax.plot(psi, np.degrees(alpha), c='b')
# ax.plot(psi, np.degrees(alpha2), c='r')
# ax.plot(psi, np.degrees(alpha3), c='g')
# plt.show()

# print('Ri:',bearing.Ri)
# print('Ro:',bearing.Ro)
# print('Rm:',bearing.dm/2)
# print('fi:',bearing.fi)
# print('fo:',bearing.fo)
# print(bearing.Display())





