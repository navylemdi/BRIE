from Material import Material
from Ball import Ball
from Raceway import Raceway
from Bearing import Bearing
import numpy as np
import matplotlib.pyplot as plt

Ceramic = Material('Ceramic')
AMS5898 = Material('AMS5898')
ball = Ball(Ceramic, 6.35/1000)
inner_raceway = Raceway(AMS5898, 'Inner', 3.397/1000, 20/1000, 12/1000, 27.959/1000)
outer_raceway = Raceway(AMS5898, 'Outer', 3.429/1000, 42/1000, 12/1000, 34.969/1000)
bearing = Bearing(outer_raceway, inner_raceway, ball, 14, 15, 30.994/1000, 0.4818/1000)
Fa=0
Fr=0
M=10
disp = bearing.solve_disp(Fa, Fr, M)
# disp_mm = bearing.variable_bar_2_variable(disp)
# print(bearing.list_variable_bar_2_variable(disp))
Qmax = bearing.Q_max(disp[0], disp[1], disp[2])
# bearing.Display()
bearing.Display_ball_load(Fa, Fr, M)
print(Qmax, 'N')
alpha0=[]
alpha=[]
alpha2=[]
alpha3=[]
psi=np.linspace(-np.pi, np.pi, bearing.Z , False)

for angle in psi:
    alpha0.append(bearing.alpha0)
    alpha.append(np.arcsin((bearing.alpha(disp[0], disp[1], disp[2], angle))))
    alpha2.append(np.arccos(bearing.alpha2(disp[0], disp[1], disp[2], angle)))
    alpha3.append(np.arctan(bearing.alpha3(disp[0], disp[1], disp[2], angle)))
# print(np.degrees(alpha))
# print(np.degrees(alpha2))
# print(np.degrees(alpha3))
# print(np.degrees(np.arcsin(bearing.alpha(0,0,0,0))))
# sum=[]
# for i in range(0, len(psi)):
#     sum.append(alpha[i]**2 + alpha2[i]**2)
# print(sum)

fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
ax.plot(psi, np.degrees(alpha0), c='k')
ax.plot(psi, np.degrees(alpha), c='b')
ax.plot(psi, np.degrees(alpha2), c='r')
ax.plot(psi, np.degrees(alpha3), c='g')
plt.show()

# print('Ri:',bearing.Ri)
# print('Ro:',bearing.Ro)
# print('Rm:',bearing.dm/2)
# print('fi:',bearing.fi)
# print('fo:',bearing.fo)
# print(bearing.Display())





