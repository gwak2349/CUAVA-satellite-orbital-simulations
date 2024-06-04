import numpy as np
import newton_method as nm
import orbitalTransforms as ot


def j2_effect(e,r_Earth,a,i,mu):
    J2 = 1.086263e-3
    Omega_dot = -(3*np.sqrt(mu)*J2*r_Earth**2)/(2*(1-e**2)**2*a**(7/2))*np.cos(i)
    
    omega_dot = Omega_dot*(5/2*np.sin(i)**2-2)/np.cos(i)
    return Omega_dot, omega_dot

def j2_on_orbit(r_Earth, i, argp, raan, Delta_t, t_f, T, mt, h, mu, e, a, t1, t_sec):
    

    Omega_dot, omega_dot = j2_effect(e, r_Earth, a, i, mu)
    M_n = 2*np.pi*t_f/T
    
    raan_i = raan + Omega_dot*Delta_t
    argp_i = argp + omega_dot*Delta_t

    

    theta, ___ = nm.calculate_anomalies(t_sec, M_n, e, mt, T)

    p_j2, q_j2, w_j2, dp_j2, dq_j2, dw_j2 = ot.elements_to_perifocal(theta, a, e, mu, h) #computing the perifocal coordinates    

    x_j2, y_j2, z_j2, dx_j2, dy_j2, dz_j2 = ot.perifocal_to_eci(p_j2, q_j2, w_j2, dp_j2, dq_j2, dw_j2, i, raan_i, argp_i) #converting perifocal to eci coordinates
    return raan_i, argp_i, x_j2, y_j2, z_j2 