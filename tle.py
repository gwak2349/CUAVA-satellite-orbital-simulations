import numpy as np
import re
import topocentric_horizon as th
import math
import j2
import ecef
import newton_method as nm
import orbitalTransforms as ot
# from sgp4.api import Satrec
# from sgp4.api import jday
import matplotlib.pyplot as plt



#Update on code: every time I change a time input into the code, the azimuth and altitude would not change in an ordered manner (i.e. by changing it by 1 second, it would go from 60 degrees to 54 degrees to 78 degrees, for example). Turns out, I had coded the sidereal time incorrectly and forgot to convert the inputs of the sinusoidal functions to radians. The code is much better now, changes in an expected manner, but the azimuth and elevation are still off, though elevation is off by about 9 degrees, basically azimuth is more of a concern since it is about 110 degrees apart.

def read_tle(filename):
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
    
    tle_elements = []

    catalogue_number = line_1[1]
    international_designator = line_1[2]
    epoch_year = float(line_1[3])
    mean_motion_dot = float(line_1[4])
    mean_motion_double_dot = line_1[5]
    B_star = line_1[6]
    ephemeis_type = line_1[7] 
    
    element_set_number = line_1[8]
    n = list(element_set_number)
    line_1_checksum = n[-1]

    tle_elements.append(catalogue_number)
    tle_elements.append(international_designator)
    tle_elements.append(epoch_year)
    tle_elements.append(mean_motion_dot)
    tle_elements.append(mean_motion_double_dot)
    tle_elements.append(B_star)
    tle_elements.append(ephemeis_type)
    tle_elements.append(element_set_number)
    tle_elements.append(line_1_checksum)
    


    inclination = float(line_2[2])*np.pi/180 #inclination
    # print(inclination*180/np.pi)
    raan = float(line_2[3])*np.pi/180 #right ascension
    eccentricity = float(line_2[4]) * 10**(-7) #eccentricity
    argument_of_perigee = float(line_2[5])*np.pi/180 #argument of perigee
    mean_anomaly = float(line_2[6])*np.pi/180 #mean anomaly
    
    y = line_2[7]
    k = list(y)
    k.remove(".")
    m = k[0:10]
    l = ''.join(map(str,m)) 
    mean_motion = float(l)*10**(-8)
     
    # mean_motion = float(line_2[7])
    n = k[10:-1]
    l1 = ''.join(map(str,n)) 
    rev_at_epoch = float(l1)
    # rev_at_epoch = line_2[8]

    # r = list(rev_at_epoch)
    line_2_checksum = int(y[-1])

    tle_elements.append(inclination)
    tle_elements.append(raan)
    tle_elements.append(eccentricity)
    tle_elements.append(argument_of_perigee)
    tle_elements.append(mean_anomaly)
    tle_elements.append(mean_motion)
    tle_elements.append(rev_at_epoch)
    # tle_elements.append(line_2_checksum)
    return tle_elements


latitude = -33.91667222
longitude = 151.033325



R = 6371e3
altitude = 391e3
position = R+altitude

filename = 'iss_tle.txt'

tle_elements = read_tle(filename)

# print(tle_elements)
R_e = 6778e3
epoch_time = tle_elements[2]
year = 2024
month = 6
day = 4
UTC = 11+8/60+0/3600

G = 6.67*10**(-11) #Gravitational constant
M = 5.97*10**24 #mass of the earth
mu = G*M # standard gravitational parameter

T = np.sqrt(R_e**3*4*np.pi**2/(G*M)) #period 
i = float(tle_elements[9])
raan = float(tle_elements[10])
e = float(tle_elements[11])
argp = float(tle_elements[12])

mt = 2*np.pi/T #initial mean motion
t0 = UTC*3600
# t1 = UTC*3600+1

ma = mt*t0 #initial mean anomaly

t_sec = np.linspace(t0,t0,1)
a = position
v = np.sqrt(mu/a)
h = a*v

raan_i, argp_i, x_j, y_j, z_j = j2.j2_on_orbit(R, i, argp, raan, 0, t0, T,mt,h,mu,e,a,t0, t_sec)
print(x_j, y_j, z_j)

# s = '1 25544U 98067A   24148.06162966  .00017600  00000-0  30994-3 0  9997'
# t = '2 25544  51.6402  60.0268 0005600 240.0999 134.8837 15.50510323455278'
# satellite = Satrec.twoline2rv(s,t)

# jd, fr = jday(2024, 5, 30, 20, 23, 25)
# e, r, v = satellite.sgp4(jd,fr) 

# x_j = [r[0]]
# y_j = [r[1]]
# z_j = [r[2]]
# print(x_j, y_j, z_j)

x_ecef, y_ecef, z_ecef, r, alpha, delta = ecef.eci_to_ecef(t0, t_sec, x_j, y_j, z_j)    

# print(alpha)
# print(delta)





# alpha = np.linspace(-90*np.pi/180,90*np.pi/180,2000)
# delta = np.linspace(-90*np.pi/180,90*np.pi/180,2000)


sidereal_time = th.calculate_sidereal_time(year, month, day, UTC, longitude)
azimuth_list = []
elevation_list = []

# for k in range(len(alpha)):
#     alpha[k] = alpha[k]*180/np.pi
#     delta[k] = delta[k]*180/np.pi

azimuth, elevation = th.geocentric_to_topocentric(alpha, delta, R_e, sidereal_time, latitude, position)
# for i in range(len(alpha)):
#     azimuth, elevation = th.geocentric_to_topocentric(alpha[i], delta[i], R_e, sidereal_time, latitude, position)
#     azimuth_list.append(azimuth)
#     elevation_list.append(elevation)

# print(azimuth*180/np.pi, elevation*180/np.pi)

print(azimuth*180/np.pi)
print(elevation*180/np.pi)
# print(azimuth_list)
# print(elevation_list)

t1 = 0

t3 = t1+86400*1000

t_sec1 = np.linspace(t1, t3, 2000)

# plt.plot(t_sec1, elevation_list)
# plt.show()