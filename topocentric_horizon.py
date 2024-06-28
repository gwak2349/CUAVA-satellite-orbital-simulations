import numpy as np
import math

def calculate_sidereal_time(y, m, d, UT, longitude):
    J0 = 367*y-math.floor(7*(y+math.floor((m+9)/12))/4)+math.floor(275*m/9)+d+1721013.5
    
    T0 = (J0-2451545)/36525
    theta_G_0 = 100.4606184+36000.77004*T0 + 0.000387933*T0**2 - 2.583*(1e-8)*T0**3
    
    if theta_G_0 > 360:
        t1 = math.floor(theta_G_0/360)
        theta = theta_G_0-t1*360
    else:
        theta = theta_G_0
    
    theta_G = theta + 360.98564724*UT/24
    sidereal_time = theta_G + longitude
    if sidereal_time > 360:
        t2 = math.floor(sidereal_time/360)
        sidereal_time = sidereal_time - t2*360
    else:
        sidereal_time = sidereal_time
   
    return sidereal_time

def geocentric_to_topocentric(alpha, delta, R_e, sidereal_time, latitude, position):
    
    r = [position*np.cos(delta)*np.cos(alpha),position*np.cos(delta)*np.sin(alpha),position*np.sin(delta)] #geocentric equatorial position vector of satellite
    R = [R_e*np.cos(latitude*np.pi/180)*np.cos(sidereal_time*np.pi/180), R_e*np.cos(latitude*np.pi/180)*np.sin(sidereal_time*np.pi/180), R_e*np.sin(latitude*np.pi/180)] #position vector of observer
    # print(R)

    rho_X = []
    for i in range(len(r)):
        c = r[i]-R[i]
        rho_X.append(c)


    # rho_X = r-R #position vector of satellite with respect to the observer in geocentric coordinates
    # print(rho_X)
    rho_x = np.array(np.matmul(geocentric_to_topocentric_matrix(latitude,sidereal_time),rho_X)) #position vector in topocentric coordinates
    # print(rho_x)
    rho_x_hat = 1/np.linalg.norm(rho_x)*rho_x
    # print(rho_x_hat)
    elevation = np.arcsin(rho_x_hat[-1])
   
    # azimuth = np.pi + np.arctan(rho_x_hat[0]/rho_x_hat[1])
    A1 = np.arcsin(rho_x_hat[0]/np.cos(elevation))
    A2 = np.arccos(rho_x_hat[1]/np.cos(elevation))

    print(A1*180/np.pi)
    # print(A2*180/np.pi)
    # if A1 > 0 and A2 > 0:
    #     azimuth = A1
    # elif A1 > 0 and A2 < 0:
    #     azimuth = A1
    # elif A1 < 0 and A2 < 0:
    #     azimuth = -A1
    # elif A1 < 0 and A2 > 0:
    #     azimuth = A2

    if A1 > 0 and A2 > 0:
        azimuth = abs(A1)
    elif A1 > 0 and A2 < 0:
        azimuth = np.pi-abs(A1)
    elif A1 < 0 and A2 < 0:
        azimuth = np.pi+abs(A1)
    elif A1 < 0 and A2 > 0:
        azimuth = 2*np.pi-abs(A1)


    return azimuth, elevation


def geocentric_to_topocentric_matrix(latitude, sidereal_time):
    
    r1 = -np.sin(sidereal_time*np.pi/180)
    r2 = np.cos(sidereal_time*np.pi/180)
    r3 = 0

    r4 = -np.sin(latitude*np.pi/180)*np.cos(sidereal_time*np.pi/180)
    r5 = -np.sin(latitude*np.pi/180)*np.sin(sidereal_time*np.pi/180)
    r6 = np.cos(latitude*np.pi/180)

    r7 = np.cos(sidereal_time*np.pi/180)*np.cos(latitude*np.pi/180)
    r8 = np.cos(latitude*np.pi/180)*np.sin(sidereal_time*np.pi/180)
    r9 = np.sin(latitude*np.pi/180)

    
    g_to_t = np.array([[r1,r2,r3],[r4,r5,r6],[r7,r8,r9]]) #geocentric to topocentric coordinates

    return g_to_t


#function below utilised concepts from https://www.geeksforgeeks.org/how-to-find-the-index-for-a-given-item-in-a-python-list/

def time_above_horizon(t,elevation):
    first_match = 10*180/np.pi
    second_match = 170*180/np.pi
    j = 0
    for i in range(len(elevation)):
        if elevation[i] > first_match:
            j = i
            break

    k = 0
    for p in range(len(elevation)):
        if elevation[p] > second_match:
            k = p
            break

    t1 = t[j]
    t2 = t[k]
    time = t2-t1
    
    return time