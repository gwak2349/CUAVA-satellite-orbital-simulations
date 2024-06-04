import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from mpl_toolkits import mplot3d
import orbitalTransforms as ot
import newton_method as nm
import re
import imageio
import geopandas as gpd
import ecef 
import j2
import julian
import read_tle as rt
import hohmann_transfer as ht

filename = "ESA_XMM_tle.txt" #TLE text file
line_1, line_2 = rt.read_tle(filename)

G = 6.6743*10**(-11) #universal gravitational constant
M = 5.9722*10**24 #Mass of the Earth

mu = G*M #standard gravitational parameter

i = float(line_2[2])*np.pi/180 #inclination
raan = float(line_2[3])*np.pi/180 #right ascension
e = float(line_2[4]) * 10**(-7) #eccentricity

argp = float(line_2[5])*np.pi/180 #argument of perigee
ma = float(line_2[6])*np.pi/180 #mean anomaly

rev_per_day = float(line_2[7])

# print(mu)
# print(i)
# print(raan)
# print(e)
# print(argp)

## Question 1 

T = (24*3600)/rev_per_day #period of orbit

mt = 2*np.pi/T #mean motion
# mt = mean_motion/(24*3600)
a = (G*M*T**2/(4*np.pi**2))**(1/3) #semi major axis

h = np.sqrt(a*mu*(1-e**2)) #specific angular momentum
r_p = a*(1-e) #perigee
r_a = 2*a-r_p #apogee

t1 = ma/mt #time at initial mean anomaly
t_sec = np.linspace(t1, T+t1, 360) #a vector of time from the initial mean anomaly until 1 revolution

# print(T)
# print(mt)
# print(a)
# print(h)
# print(r_p)
# print(r_a)

# position = []
# velocity = []
mean_anomaly, theta = nm.calculate_anomalies(t_sec, ma, e, mt, T)

# for ta in theta:
            
#     #True anomaly     

#     r = h**2/(mu*(1+(e*np.cos(ta))))
#     v = mu/h*(e+np.cos(ta))
    
#     position.append(r)
#     velocity.append(v)

p, q, w, dp, dq, dw = ot.elements_to_perifocal(theta, a, e, mu, h) #computing the perifocal coordinates
# print(e)
p_to_e = ot.perifocal_to_eci_matrix(i, raan, argp) #obtaining the rotation matrix to convert perifocal to eci coordinates

x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp) #converting perifocal to eci coordinates


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

r_Earth = 6.371*10**6

# #code below was adjusted from https://www.tutorialspoint.com/plotting-points-on-the-surface-of-a-sphere-in-python-s-matplotlib
u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
x_Earth = r_Earth*np.cos(u) * np.sin(v)
y_Earth = r_Earth*np.sin(u) * np.sin(v)
z_Earth = r_Earth*np.cos(v)

# ax.plot(x,y,z, color="green") #plots the orbit in eci coordinates
# ax.plot_surface(x1, y1, z1, color="blue") #plots the spherical Earth
# # plt.title("Plot of satellite's orbit in eci coordinates")
# ax.set_xlabel('x-coordinate (m) 1e6',fontsize=25,labelpad=25)
# ax.set_ylabel('y-coordinate (m) 1e6',fontsize=25,labelpad=25)
# ax.set_zlabel('z-coordinate (m) 1e6',fontsize=25,labelpad=25)
# ax.set_aspect('equal', adjustable='box')


# ax.tick_params(axis='x', labelsize=25)
# ax.tick_params(axis='y', labelsize=25)
# ax.tick_params(axis='z', labelsize=25)


# # ax.xaxis.label.set_size(50)
# # ax.yaxis.label.set_size(50)
# # ax.zaxis.label.set_size(50)


# # ax.xticks(size = 15)
# # ax.yticks(size = 15)
# # ax.zticks(size = 15)
# plt.show()

t_epoch = line_1[3]


#https://stackoverflow.com/questions/2579535/convert-dd-decimal-degrees-to-dms-degrees-minutes-seconds-in-python

sidereal_time = julian.sidereal_time(t_epoch)

x_ecef, y_ecef, z_ecef, r, alpha, delta = ecef.eci_to_ecef(sidereal_time,t_sec,x,y,z)

# plot of the map was generated using the code from https://towardsdatascience.com/the-easiest-way-to-plot-data-from-pandas-on-a-world-map-1a62962a27f3
countries = gpd.read_file(
               gpd.datasets.get_path("naturalearth_lowres"))
countries.plot(color="lightgrey")
plt.plot(alpha,delta,'.')
# plt.xlabel("Right ascension of the ascending node (degrees) of XMM-Newton satellite", fontsize = 25, labelpad = 10)
# plt.ylabel("Declination of XMM-Newton satellite (degrees)",fontsize = 25, labelpad = 10)
# plt.xticks(size = 25)
# plt.yticks(size = 25)

# plt.grid()
plt.show()

deltat_1 = 60*60*24*2
# deltat_n = 60*60*24*4
deltat_n = 60*60*24*2*7*52*12
delta_t = np.linspace(deltat_1,deltat_n,5)

time_1 = t1+deltat_1
time_2 = t1+deltat_n

time = np.linspace(time_1, time_2, 5)

Omega_dot, omega_dot = j2.j2_effect(e,r_Earth,a,i,mu)
for t in range(len(time)):
    x_j2, y_j2, z_j2, x_ecef_j2, y_ecef_j2, z_ecef_j2, r_j2, alpha_j2, delta_j2, mean_anomaly = j2.j2_on_ecef(Omega_dot, omega_dot, i, argp, raan, delta_t[t],time[t],T,mt,h,mu,e,a)


#     ax.plot(x_j2, y_j2, z_j2, linestyle='solid') 
# ax.plot_surface(x_Earth, y_Earth, z_Earth, color="blue") #plots the spherical Earth
# ax.set_xlabel('x-coordinate (m) 1e6',fontsize=25, labelpad = 25)
# ax.set_ylabel('y-coordinate (m) 1e6',fontsize=25, labelpad = 25)
# ax.set_zlabel('z-coordinate (m) 1e6',fontsize=25, labelpad = 25)

# ax.tick_params(axis='x', labelsize=15)
# ax.tick_params(axis='y', labelsize=15)
# ax.tick_params(axis='z', labelsize=15)
# ax.legend(['1 period', '3 years later', '6 years later', '9 years later', '12 years later'],fontsize=15,loc='upper left')
# ax.set_aspect('equal', adjustable='box')
# plt.show()
    
# countries = gpd.read_file(
#             gpd.datasets.get_path("naturalearth_lowres"))    
# countries.plot(color="lightgrey")
for t in range(len(time)):
    
    x_j2, y_j2, z_j2, x_ecef_j2, y_ecef_j2, z_ecef_j2, r_j2, alpha_j2, delta_j2, mean_anomaly = j2.j2_on_ecef(Omega_dot, omega_dot, i, argp, raan, delta_t[t],time[t],T,mt,h,mu,e,a)
    
#     plt.plot(alpha_j2,delta_j2,'.', markersize=1)
# plt.show()


## Question 2


i1 = 0
r1 = 500*1e3
rp_1 = r1
ra_1 = r1
a1 = (rp_1+ra_1)/2
v1 = np.sqrt(mu/r1)
h1= r1*v1
e1 = 0
T1 = np.sqrt(r1**(3/2)*4*np.pi**2/mu)
raan1 = 0
argp1 = 0
t_sec1 = np.linspace(t1, T1+t1, 360)
mt1 = 2*np.pi/T1
ma1 = mt1*t_sec[0]
mean_anomaly1, theta1 = nm.calculate_anomalies(t_sec1, ma1, e1, mt1, T1)
p1, q1, w1, dp1, dq1, dw1 = ot.elements_to_perifocal(theta1, r1, e1, mu, h1) #computing the perifocal coordinates

x1, y1, z1, dx1, dy1, dz1 = ot.perifocal_to_eci(p1, q1, w1, dp1, dq1, dw1, i1, raan1, argp1) #converting perifocal to eci coordinates
# delta_vA, delta_vB, t_manoeuvre = ht.find_transfer_orbit(rp_1, ra_1, e1, h1, a1, a3, rp_3, ra_3, e3, h3, mu, T1, T3, mt, t_sec)
# plt.plot(p,q)
# plt.plot(p1,q1)
# plt.axhline(y=0, color='r', linestyle='-')
# plt.axis('square')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# ax.plot(x,y,z)
# ax.plot(x1,y1,z1)
# plt.show()

# p_to_e = ot.perifocal_to_eci_matrix(i1, raan1, argp1) #obtaining the rotation matrix to convert perifocal to eci coordinates
# a3 = a
# rp_3 = r_p
# ra_3 = r_a
# h3 = h
# T3 = T
# e3 = 3

# delta_vA, delta_vB, t_manoeuvre, i2, raan2, argp2 = ht.find_transfer_orbit(rp_1, ra_1, e1, h1, a1, a3, rp_3, ra_3, e3, h3, mu, T1, T3, mt, t_sec, ma)




# print(p1)
# print(q1)


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.plot(x,y,z)
# ax.plot(x1,y1,z1)
# ax.set_aspect('equal', adjustable='box')
# ax.set_box_aspect([1,1,1])
# ax = plt.gca()
# ax.set_aspect('equal', adjustable='box')
# plt.show()