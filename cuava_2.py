import numpy as np 
import newton_method as nm
import orbitalTransforms as ot
import matplotlib.pyplot as plt
import j2
import ecef
import geopandas as gpd
import tle
import eci_plot as ep
import groud_trace_plot as gtp

altitude = 510e3 #altitude
e = 0 #eccentricity
R = 6371e3 #radius of the earth
G = 6.67*10**(-11) #Gravitational constant
a = altitude+R #distance from the centre of the Earth to the satellite, i.e. semi-major axis
M = 5.97*10**24 #mass of the earth
mu = G*M # standard gravitational parameter

T = np.sqrt(a**3*4*np.pi**2/(G*M)) #period of the satellite's orbit

Omega_dot = 1.991*10**(-7) #rate of change of RAAN
J2 = 1.08263*10**(-3) #J2 coefficient
i = np.arccos((Omega_dot*-2*(1-e**2)**2*a**(7/2))/(3*np.sqrt(mu)*J2*R**2)) #inclination for sun synchronous orbit
raan = 337.5*np.pi/180 #right ascension of the ascending node
argp = 0 #argument of perigee

mt = 2*np.pi/T #mean motion of satellite
ma = 0 #initial mean anomaly
t1 = ma/mt #initial time

t3 = 4*T #4 orbits
t_sec = np.linspace(t1,t3,5000) #time array from t=0 to 4 orbits
v = np.sqrt(mu/a) #orbital speed of satellite
h = a*v #angular momentum of satellite

mean_anomaly, theta = nm.calculate_anomalies(t_sec, ma, e, mt, T) #calculating mean anomaly and true anomaly vector of the satellite
p, q, w, dp, dq, dw = ot.elements_to_perifocal(theta, a, e, mu, h) #computing the perifocal coordinates of the satellite

x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp) #converting perifocal coordinates to eci coordinates

r_Earth = 6.371*10**6 #radius of the Earth

## Calculating right ascension and declination of satellite to create a ground trace plot
raan_i, argp_i, x_j, y_j, z_j = j2.j2_on_orbit(r_Earth, i, argp, raan, 86400, t3, T,mt,h,mu,e,a,t1, t_sec)
x_ecef, y_ecef, z_ecef, r, alpha, delta = ecef.eci_to_ecef(0, t_sec, x_j, y_j, z_j)    

## Ground trace plot
gtp.ground_trace_plot(alpha, delta)







