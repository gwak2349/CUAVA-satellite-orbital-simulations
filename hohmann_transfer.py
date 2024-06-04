import numpy as np
import orbitalTransforms as ot
import newton_method as nm

def find_transfer_orbit(rp_1, ra_1, e1, h1, a1, a3, rp_3, ra_3, e3, h3, mu, T1, T3, mt, t_sec, ma, i2, raan2, argp2):


    rp_2 = rp_1 #perigee of transfer orbit
    ra_2 = ra_3 #apogee of transfer orbit
    a2 = (ra_2+rp_2)/2 #semi-major axis of transfer orbit
    e2 = (ra_2-rp_2)/(ra_2+rp_2) #eccentricity of transfer orbit
    h2 = np.sqrt(2*mu*ra_2*rp_2/(ra_2+rp_2)) #specific angular momentum of transfer orbit
    T2 = np.sqrt((a2**3*4*np.pi**2)/mu) #period of transfer orbit

    v_A1 = h1/rp_1 #velocity at the perigee of the initial orbit
    v_A2 = h2/rp_1 #velocity at the perigee of the transfer orbit
    v_B1 = h3/ra_2 #velocity at the apogee of the transfer orbit
    v_B2 = h2/ra_2 #velocity at the apogee of the final orbit

    delta_vA = v_A2 - v_A1 #change in velocity to go from initial orbit to transfer orbit
    delta_vB = v_B2- v_B1 #change in velocity to go from transfer orbit to final orbit

    t_manoeuvre = 0.5*T2

    t_sec2 = np.linspace(t_manoeuvre, t_manoeuvre+T2, 1000) #time vector starting from the perigee to the apogee of the transfer orbit
    ma2 = np.pi #mean anomaly at initial time of transfer orbit
    mt2 = 2*np.pi/T2 #mean motion of the transfer orbit
    theta2, ea2, ma_vec2 = ot.compute_anomalies(t_sec2, ma2, e2, mt2) #computing the mean, eccentric and true anomalies at every timestep in the t_sec2 vector
    p2, q2, w2, dp2, dq2, dw2 = ot.elements_to_perifocal(theta2, a2, e2, mu, h2) #calculating the perifocal coordinates of the orbit from perigee to apogee of tranfer orbit

    print(len(p2),len(q2))


    # x2, y2, z2, dx2, dy2, dz2 = ot.perifocal_to_eci(p2, q2, w2, dp2, dq2, dw2, i2, raan2, argp2) #converting perifocal to eci coordinates
    # return  delta_vA, delta_vB, t_manoeuvre, x2, y2, z2, dx2, dy2, dz2
    return  delta_vA, delta_vB, t_manoeuvre

# def calculate_orbital_elements(position, velocity, mu, e2):
    h_vector = []
    for i in range(len(p2)):
    
        sam = np.matmul(position,velocity) #specific angular momentum
        h_vector.append(sam)

    print(h_vector)
    h_x = h_vector[0]
    h_y = h_vector[1]
    h_z = h_vector[2]
    h = np.sqrt(h_x**2+h_y**2+h_z**2)
    i = np.arccos(h_z/h)

    K = [0,0,1]
    N_vector = np.array(np.matmul(K,h_vector))
    N_x = N_vector[0]
    N_y = N_vector[1]
    N_z = N_vector[2]
    N = np.sqrt(N_x**2+N_y**2+N_z**2)
    if N_y < 0:
        raan = 2*np.pi-np.arccos(N_x/N)
    else:
        raan = np.arccos(N_x/N)

    e_vector = 1/mu*(np.linalg.norm(velocity)-mu/np.linalg.norm(position)*position-np.linalg.norm(position)*v_r*velocity)
    e_x = e_vector[0]
    e_y = e_vector[1]
    e_z = e_vector[2]
    
    if e_z < 0:
        argp = 2*np.pi-np.arccos(np.dot(N_vector, e_vector)/(N*e2))
    else:
        argp = np.arccos(np.dot(N_vector, e_vector)/(N*e2))
    return i, raan, argp






# mean_anomaly, theta = nm.calculate_anomalies(t_sec, ma, e1, mt, T1)


# true_anomaly = theta[0] #true anomaly at the initial epoch time
# A = np.sqrt((1-e1)/(1+e1))
# E = 2*np.arctan(A*np.tan(true_anomaly/2)) #eccentric anomaly at the initial epoch time
# M_e = E-e1*np.sin(E) #mean anomaly at the initial epoch time
# t_start = (T1*M_e/(2*np.pi)) #initial time
# t_total = t_manoeuvre+(T1-t_start)

# t_sec1 = np.linspace(t_start, T1 ,150) #vector that starts from the initial epoch time until the time at the perigee

# ta, ea, ma_vec = ot.compute_anomalies(t_sec1, ma, e1, mt) #computing the mean, eccentric and true anomalies at every timestep in the t_sec1 vector

# p1, q1, w1, dp1, dq1, dw1 = ot.elements_to_perifocal(ta, a1, e1, mu, h1) #calculating the perifocal coordinates of the orbit from start time to perigee




# t_sec3 = np.linspace(-t_manoeuvre,-T3, 180) #time vector starting from the perigee to the apogee of the transfer orbit
# ma3 = ma_vec[-1] + ma #mean anomaly at initial time of transfer orbit
# mt3 = 2*np.pi/T3 #mean motion of the transfer orbit
# ta3, ea3, ma_vec3 = ot.compute_anomalies(t_sec3, ma3, e3, mt3) #computing the mean, eccentric and true anomalies at every timestep in the t_sec3 vector
# p3, q3, w3, dp3, dq3, dw3 = ot.elements_to_perifocal(ta3, a3, e3, mu, h3) #calculating the perifocal coordinates of the orbit from perigee to apogee of tranfer orbit


# t_sec2 = np.linspace(t_manoeuvre, t_manoeuvre+T2, 1000) #time vector showing one period of the final orbit
# ma2 = np.pi #initial mean anomaly of final orbit
# mt3 = 2*np.pi/T3 #mean motion of final orbit
# ta2, ea2, ma_vec2 = ot.compute_anomalies(t_sec2, ma2, e3, mt3) #computing the mean, eccentric and true anomalies at every timestep in the t_sec2 vector
# p2, q2, w2, dp2, dq2, dw2 = ot.elements_to_perifocal(ta2, a2, e2, mu, h2) #calculating the perifocal coordinates of final orbit



