import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize
import plotly.graph_objects as go
import time

#Ceramic
E1 = 315_000e6 #Material 1 Young's modulus [Pa]
D11 = 6.35e-3 #Ball 1 Diameter [mm]
D12 = D11
Nu1 = 0.26 #Poisson's ratio 1

#AMS5898
E2 = 208_000e6 #Material 2 Young's modulus [Pa]
D21 = 25e-3 #raceway 1 diameter [mm]
D22 = -10e-3 #raceway 2 diameter [mm]
Nu2 = 0.3 #Poisson's ratio 2

p=100 #Force [N]
M_TO_MM = 1000
MM_TO_M = 1/1000
N_TO_KG = 1/9.81
KG_TO_N = 9.81

def convert_time(time):
    hour=time//3600
    if time//3600==0:
        minutes=time//60
        if minutes//60!=0:
            second=time
        else:
            second = time-minutes*60
    else:
        minutes = time//60 - hour*60
        second = time-minutes*60-hour*60*60
    print("Elapsed time:", hour, 'h', minutes,'m', np.round(second), 's')

def theta(nu):
    return nu/(1-2*nu)

def k(E, nu):
    return E/(1+nu)/2

def contraction_theta(k, theta):
    return (2*(1+theta))/(k*(1+2*theta))

contraction_theta1 = contraction_theta(k(E1, Nu1), theta(Nu1))
contraction_theta2 = contraction_theta(k(E2, Nu2), theta(Nu2))

print('\u03D11=',contraction_theta1, '\u03D12=', contraction_theta2)

def A(D11, D12, D22, D21):
    rho11, rho12, rho22, rho21 = 2/D11, 2/D12, 2/D22, 2/D21
    return 1/4 * np.sqrt((rho11-rho12)**2+2*(rho11-rho12)*(rho21-rho22)+(rho21-rho22)**2)

def B(D11, D12, D22, D21):
    rho11, rho12, rho22, rho21 = 1/D11, 1/D12, 1/D22, 1/D21
    return 1/2 * (rho11+rho12+rho21+rho22) - 1/2*np.sqrt((rho11-rho12)**2+2*(rho11-rho12)*(rho21-rho22)+(rho21-rho22)**2)
Var_A=A(D11, D12, D22, D21)
Var_B=B(D11,D12,D22,D21)

print('A=', A(D11, D12, D22, D21), 'B=', B(D11,D12,D22,D21))

K1 = Var_A * 8* np.pi/((contraction_theta1+contraction_theta2)*3*p)
K2 = Var_B * 8* np.pi/((contraction_theta1+contraction_theta2)*3*p)

# def F(lamb):
#     integrand =lambda t:1/np.sqrt((1+t)**3*(lamb**2+t)*t)
#     result, abserr = integrate.quad(integrand, 0, np.inf, limit=200)
#     # print('Absolute error estimated', abserr)
#     return result

# def H(lam, K1, K2):
#     return F(1/lam)/(lam**3 * F(lam)) - K2/K1

# def solve_ab(K1, K2):

#     # résolution pour lambda
#     lam_solution = optimize.root_scalar(
#         H,
#         args=(K1, K2),
#         bracket=[1e-3, 1e3],   # intervalle large
#         method='brentq'
#     ).root

#     # calcul de a
#     a = (F(lam_solution)/K1)**(1/3)

#     # calcul de b
#     b = a * lam_solution

#     return a, b

# def LoadvsEllipse():
    # p=np.linspace(1, 10000, 100)
    # a_array=[]
    # b_array=[]
    # for load in p:
    #     K1 = Var_A * 16* np.pi/((contraction_theta1+contraction_theta2)*3*load)
    #     K2 = Var_B * 16* np.pi/((contraction_theta1+contraction_theta2)*3*load)
    #     a_array.append(solve_ab(K1,K2)[0]*M_TO_MM)
    #     b_array.append(solve_ab(K1,K2)[1]*M_TO_MM)
    # fig=plt.figure()
    # plt.plot(p, a_array, label='Major semi axis a')
    # plt.plot(p, b_array, label='Minor semi axis b')
    
    # plt.xlabel('Load [N]')
    # plt.ylabel('Length [mm]')
    # plt.title(f'Ellipse axes length versus load \n for a {D11*1000}mm ceramic ball \n and a ({D21*1000}, {D22*1000})mm inner bearing AMS5898 raceway')
    # plt.grid()
    # # plt.show()

def F2(k):
    integrand =lambda t: 1/np.sqrt((1+(k*t)**2)**3 * (1+t**2))
    result, abserr = integrate.quad(integrand, 0, np.inf, limit=200)
    # print('Absolute error estimated', abserr)
    return result

def H2(k, K1, K2):
    return (F2(k) * k**2 / F2(1/k)) - K1/K2

def solve_ab2(K1, K2):

    # résolution pour k
    k_solution = optimize.root_scalar(
        H2,
        args=(K1, K2),
        bracket=[1e-3, 1e3],   # intervalle large
        method='brentq'
    ).root

    # calcul de a
    a = (F2(k_solution)/K1)**(1/3)

    # calcul de b
    b = a / k_solution

    return a, b


def distance_a(k, contraction_theta1, contraction_theta2, p, semis_axes_a):
    return F2(k)*3*p*(contraction_theta1+contraction_theta2)/(semis_axes_a*8*np.pi)

def LoadvsEllipse2():
    load=np.linspace(1, 5000, 50)
    a_array=[]
    b_array=[]
    distance_a_array=[]
    for p in load:
        k1 = Var_A * 8* np.pi/((contraction_theta1+contraction_theta2)*3*p)
        k2 = Var_B * 8* np.pi/((contraction_theta1+contraction_theta2)*3*p)
        a,b=solve_ab2(k1,k2)        
        k=a/b
        a_array.append(a*M_TO_MM)
        b_array.append(b*M_TO_MM)
        distance_a_array.append(distance_a(k, contraction_theta1, contraction_theta2, p, a)*M_TO_MM)
    
    plt.plot(load, a_array, label='Major semi axis a_2')
    plt.plot(load, b_array, label='Minor semi axis b_2')
    plt.plot(load, distance_a_array, label='Distance of penetration')
    plt.legend()
    plt.xlabel('Load [N]')
    plt.ylabel('Length [mm]')
    plt.title(f'Ellipse axes length versus load \n for a {D11*1000}mm ceramic ball \n and a ({D21*1000}, {D22*1000})mm inner bearing AMS5898 raceway')
    plt.grid()
    plt.show()
    
# LoadvsEllipse()

# LoadvsEllipse2()

def Hu(u, a,b,x,y,z):
    return x**2/(a**2+u) + y**2/(b**2+u) + z**2/u -1

def solve_u(a,b,x,y,z):
    if (x/a)**2 + (y/b)**2< 1:
        return 0
    else:
        # try:
        #     results = optimize.root_scalar(
        #     Hu,
        #     args=(a,b,x,y,z),
        #     bracket=[1e-6, 1e3],   # intervalle large
        #     method='brentq'
        #     ).root
        #     return results
        # except ValueError:
        #     print(a,b,x,y,z)
        A=1
        B=(a**2+b**2-x**2-y**2-z**2)
        C=(a**2*b**2-x**2*b**2-y**2*a**2-z**2*a**2-z**2*b**2)
        D=-z**2*a**2*b**2
        return max(np.roots([A,B,C,D]))
    

    
a,b = solve_ab2(K1, K2)
print('Major axes', a, "m,", a*M_TO_MM, "mm", a*M_TO_MM*M_TO_MM, "µm",)
print('Minor axes', b, "m,", b*M_TO_MM, "mm", b*M_TO_MM*M_TO_MM, "µm",)
# x,y,z=1,1,0
# u = solve_u(a,b,x,y,z)
def contour_ellipse(a,b,N):
    theta=np.linspace(0,2*np.pi, N, endpoint=True)
    r=a*b/(np.sqrt((a*np.sin(theta))**2 + (b*np.cos(theta))**2))
    x_ellipse = r*np.cos(theta)
    y_ellipse = r*np.sin(theta)
    # y_ellipse=np.sqrt((1-(x/a)**2))*b
    # print(x)
    # print((x/a)**2)
    # print(y_ellipse)
    return x_ellipse, y_ellipse

def P(p, a, b, u, x, y, z):
    integrand= lambda t: (1 - x**2/(a**2 + t) - y**2/(b**2 + t) - z**2/t) / np.sqrt((a**2+t)*(b**2+t)*t)
    result, abserr = integrate.quad(integrand, u, np.inf, limit=200)
    return result * 3*p/(16*np.pi)

def Zeta(x,y,z, contraction_theta, load, P):
    integrand= lambda t: 1/(t*np.sqrt((a**2+t)*(b**2+t)*t))
    integration, abserr = integrate.quad(integrand, u, np.inf, limit=200)
    return 6*load*z**2*integration/(16*np.pi) + contraction_theta*P

# print(P(p,a,b,u,x,y,z))
N_resolution = 10
if N_resolution%2 == 0 :
    N_resolution=N_resolution+1

x=np.linspace(-max(a,b), max(a,b), N_resolution)
y=np.linspace(-max(a,b), max(a,b), N_resolution)
z=np.linspace(-max(a,b), max(a,b), N_resolution)
x_limited=np.linspace(0, max(a,b), N_resolution)
y_limited=np.linspace(0, max(a,b), N_resolution)
z_limited=np.linspace(0, max(a,b), N_resolution)

compt=0
# print('Estimated time of resolution', 2*N_resolution**3 , 's')
s=time.time()
P_array=np.empty((N_resolution,N_resolution,N_resolution))
Zeta1_array=np.empty((N_resolution,N_resolution,N_resolution))

for id_x, x_position in enumerate(x):
    for id_y, y_position in enumerate(y):
        for id_z, z_position in enumerate(z):
            # print(x_position,y_position,z_position)
            u = solve_u(a,b,x_position,y_position,z_position)
            # print(u)
            compt+=1
            print (np.round(compt/N_resolution**3*100,2), '%', end="\r")
            P_array[id_x, id_y, id_z] = P(p,a,b,u,x_position,y_position,z_position)
            Zeta1_array[id_x, id_y, id_z] = Zeta(x_position, y_position, z_position, contraction_theta1, p, P_array[id_x, id_y, id_z])



print('Resolution : ', N_resolution,' nodes')
convert_time(time.time()-s) #Affichage de la durée de calcul des boucles for

u = solve_u(a,b,0,0,0)
print('Z displacement:', Zeta(0,0,0, contraction_theta1, p, P(p,a,b,u,1,1,1))*1000, 'mm')
X,Y,Z = np.meshgrid(x,y,z, indexing='ij')
Y_surface, Z_surface = np.meshgrid(y,z, indexing='ij')# pour slice_x
X_surface, Z2_surface = np.meshgrid(x,z, indexing='ij')# pour slice_y
X2_surface, Y2_surface = np.meshgrid(x,y, indexing='ij')# pour slice_z
z1 = 1/D11*X_surface**2 + 1/D12*Y_surface**2
z2 = 1/D22*X_surface**2 + 1/D21*Y_surface**2

x_slice = 0
index_x_slice = np.argwhere(x==x_slice)
y_slice = 0
index_y_slice = np.argwhere(y==y_slice)
z_slice = 0
index_z_slice = np.argwhere(z==z_slice)

P_array_x_slice = np.reshape(P_array[index_x_slice,:,:],(N_resolution,N_resolution))
P_array_y_slice = np.reshape(P_array[:,index_y_slice,:],(N_resolution,N_resolution))
P_array_z_slice = np.reshape(P_array[:,:,index_z_slice],(N_resolution,N_resolution))
Zeta1_array_x_slice = np.reshape(Zeta1_array[index_x_slice,:,:],(N_resolution,N_resolution))
Zeta1_array_y_slice = np.reshape(Zeta1_array[:,index_y_slice,:],(N_resolution,N_resolution))
Zeta1_array_z_slice = np.reshape(Zeta1_array[:,:,index_z_slice],(N_resolution,N_resolution))
# print(X_surface,X_surface.shape)
# print(X2_surface.shape)
# print(Y_surface.shape)
# print(Y2_surface.shape)
# print(Z_surface.shape)
# print(P_array_x_slice, P_array_x_slice.shape)
# print(P_array_y_slice, P_array_y_slice.shape)
# print(P_array_z_slice, P_array_z_slice.shape)


fig, axes = plt.subplots(2,2)

pcm0 = axes[0][0].pcolormesh(Y2_surface, X2_surface, P_array_z_slice)
axes[0][0].set_title('Z_slice')
axes[0][0].set_xlabel('Y [mm]')
axes[0][0].set_ylabel('X [mm]')
axes[0][0].plot(contour_ellipse(a, b, N_resolution**3)[1],contour_ellipse(a, b, N_resolution**3)[0], c='b')
# axes[0][0].plot(-contour_ellipse(a,b,x_ellipse),x_ellipse, c='b')



pcm1 = axes[0][1].pcolormesh(Z2_surface, X_surface, P_array_y_slice)
axes[0][1].set_title('Y_slice')
axes[0][1].set_xlabel('Z [mm]')
axes[0][1].set_ylabel('X [mm]')

pcm2 = axes[1][0].pcolormesh(Y_surface, Z_surface, P_array_x_slice)
axes[1][0].set_title('X_slice')
axes[1][0].set_xlabel('Y [mm]')
axes[1][0].set_ylabel('Z [mm]')


pcm3 = axes[1][1].pcolormesh(Y_surface, Z_surface, Zeta1_array_x_slice)
axes[1][1].set_title('X_slice')
axes[1][1].set_xlabel('Y [mm]')
axes[1][1].set_ylabel('Z [mm]')
fig.colorbar(pcm0)
fig.colorbar(pcm1)
fig.colorbar(pcm2)
fig.colorbar(pcm3)
plt.show()

# fig = go.Figure()
# fig.add_trace(go.Volume(
#     x=X.flatten(),
#     y=Y.flatten(),
#     z=Z.flatten(),
#     value=P_array.flatten(),
#     isomin=np.amin(P_array),
#     isomax=np.amax(P_array),
#     opacity=0.1, # needs to be small to see through all surfaces
#     surface_count=21, # needs to be a large number for good volume rendering
#     ))
# fig.add_trace(go.Surface(
#     x=X_surface,
#     y=Y_surface,
#     z=z1,
#     surfacecolor=np.ones_like(z1),  # champ constant
#     colorscale=[[0, 'black'], [1, 'black']],
#     ))

# fig.add_trace(go.Surface(
#     x=X_surface,
#     y=Y_surface,
#     z=z2,
#     surfacecolor=np.ones_like(z2),  # champ constant
#     colorscale=[[0, 'grey'], [1, 'grey']],))
# fig.update_layout(
#     scene=dict(
#         zaxis=dict(range=[-1,1])  # limite visuelle
#     )
# )
# fig.show()
