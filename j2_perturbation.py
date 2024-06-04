import numpy as np
import math
import newton_method as nm
import orbitalTransforms as ot
import ecef 
import re
import julian
import read_tle as rt
import matplotlib.pyplot as plt


def j2_effect(e,r_Earth,a,i,mu, argp, raan, Delta_t,t_f,T,mt,h):
    J2 = 1.086263e-3
    Omega_dot = -(3*np.sqrt(mu)*J2*r_Earth**2)/(2*(1-e**2)**2*a**(7/2))*np.cos(i)
    
    omega_dot = Omega_dot*(5/2*np.sin(i)**2-2)/np.cos(i)

    t1 = t_f - Delta_t
    t_sec = np.linspace(t1, t_f, 1000)
    M_n = 2*np.pi*t1/T
    theta, mean_anomaly = nm.calculate_anomalies(t_sec, M_n, e, mt, T)
    p_j2, q_j2, w_j2, dp_j2, dq_j2, dw_j2 = ot.elements_to_perifocal(theta, a, e, mu, h) #computing the perifocal coordinates
    # print(e)
    # p_to_e = ot.perifocal_to_eci_matrix(i, Omega_n, w_n) #obtaining the rotation matrix to convert perifocal to eci coordinates
    
    filename = "ESA_XMM_tle.txt" #TLE text file
    line_1,line_2 = rt.read_tle(filename)
    
    t_epoch = line_1[3]
    
    sidereal_time = julian.sidereal_time(t_epoch)

    Omega_n = raan + Omega_dot*Delta_t
    
    w_n = argp + omega_dot*Delta_t

    x_j2, y_j2, z_j2, dx_j2, dy_j2, dz_j2 = ot.perifocal_to_eci(p_j2, q_j2, w_j2, dp_j2, dq_j2, dw_j2, i, Omega_n, w_n) #converting perifocal to eci coordinates
    x_ecef_j2, y_ecef_j2, z_ecef_j2, r_j2, alpha_j2, delta_j2 = ecef.eci_to_ecef(sidereal_time,t_sec,x_j2,y_j2,z_j2)
    
    return x_j2, y_j2, z_j2, x_ecef_j2, y_ecef_j2, z_ecef_j2, r_j2, alpha_j2, delta_j2, mean_anomaly

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

mean_motion = float(line_2[7])



## Question 1 

r_Earth = 6.371*10**6
T = (24*3600)/mean_motion #period of orbit
# print(T)
mt = 2*np.pi/T #mean motion
# mt = mean_motion/(24*3600)
a = (G*M*T**2/(4*np.pi**2))**(1/3) #semi major axis

h = np.sqrt(a*mu*(1-e**2)) #specific angular momentum
t1 = ma/mt



deltat_1 = 60*60*24*2
deltat_n = 60*60*24*2*7*52*12
delta_t = np.linspace(deltat_1,deltat_n,5)



time_1 = t1+deltat_1
time_2 = t1+deltat_n

time = np.linspace(time_1, time_2, 5)



for t in range(len(time)):
    x_j2, y_j2, z_j2, x_ecef_j2, y_ecef_j2, z_ecef_j2, r_j2, alpha_j2, delta_j2, mean_anomaly = j2_effect(e,r_Earth,a,i,mu, argp, raan, delta_t[t],time[t],T,mt,h)
    
    
    plt.plot(alpha_j2,delta_j2,'.', markersize=1)
plt.show()
