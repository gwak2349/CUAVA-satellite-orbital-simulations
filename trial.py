import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from mpl_toolkits import mplot3d
import orbitalTransforms as ot


filename = "ESA_XMM_tle.txt" #TLE text file
txt = open(filename, "r") #opening the TLE text file
contents = txt.read() #reading the TLE text file
TLE = []

for x in contents.split("\n"):
    TLE.append(x.split(" "[-1])) #splitting each element in the TLE  
    
line_0 = TLE[0]
line1 = TLE[1]
line2 = TLE[2]

line_1 = [x for x in line1 if x != '']
line_2 = [x for x in line2 if x != '']

G = 6.6743*10**(-11) #universal gravitational constant
M = 5.9722*10**24 #Mass of the Earth

mu = G*M #standard gravitational parameter

i = float(line_2[2])*np.pi/180 #inclination
raan = float(line_2[3])*np.pi/180 #right ascension
e = float(line_2[4]) * 10**(-7) #eccentricity

argp = float(line_2[5])*np.pi/180 #argument of perigee
ma = float(line_2[6])*np.pi/180 #mean anomaly

mean_motion = float(line_2[7])


##Question 1

T = (24*3600)/mean_motion #period of orbit

mt = 2*np.pi/T #mean motion
# mt = mean_motion/(24*3600)
a = (G*M*T**2/(4*np.pi**2))**(1/3) #semi major axis

h = np.sqrt(a*mu*(1-e**2)) #specific angular momentum
r_p = a*(1-e) #perigee
r_a = 2*a-r_p #apogee

t1 = ma/mt #time at initial mean anomaly
t_sec = np.linspace(t1, T+t1, 360) #a vector of time from the initial mean anomaly until 1 revolution
ta, ea, ma_vec = ot.compute_anomalies(t_sec, ma, e, mt) #calculating the true, eccentric and mean anomaly

v_r, v_n = ot.compute_orbital_velocity(a, e, ta, mu) #computing the orbital velocities
# print(e)
p, q, w, dp, dq, dw = ot.elements_to_perifocal(ta, a, e, mu, h) #computing the perifocal coordinates
# print(e)
p_to_e = ot.perifocal_to_eci_matrix(i, raan, argp) #obtaining the rotation matrix to convert perifocal to eci coordinates

x, y, z, dx, dy, dz = ot.perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp) #converting perifocal to eci coordinates



position = []
velocity = []

for i in range(len(p)):
    r = np.sqrt(p[i]**2+q[i]**2)
    position.append(r) #appending each perifocal distance to the position list
  

for i in range(len(v_r)):
    v = np.sqrt(v_r[i]**2+v_n[i]**2)
    velocity.append(v) ##appending each velocity value to the velocity list


plt.plot(t_sec, position) #Plotting position vs time of the satellite 
plt.title("Position of satellite vs time")
plt.xlabel("Time (s)")
plt.ylabel("Position(m)")
plt.show()

# plt.plot(t_sec,velocity) #Plotting velocity vs time of the satellite
# plt.title("Velocity of satellite vs time")
# plt.xlabel("Time(s)")
# plt.ylabel("Velocity(m/s)")
# plt.show()


# plt.plot(p,q)
# plt.title("Plot of satellite's orbits in perifocal coordinates")
# plt.xlabel("Perifocal coordinate p (m)")
# plt.ylabel("Perifocal coordinate q (m)")

# plt.show()