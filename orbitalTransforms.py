import numpy as np
import scipy as sp


def compute_orbital_velocity(a, e, ta, mu):
    #Calculating orbital velocities
    v_r = []
    v_n = []
    for i in ta:
        vr = float(np.sqrt(mu/a*(1+e*np.cos(i)))*e*np.sin(i))  # Radial velocity
        vn = float(np.sqrt((mu*(1+e*np.cos(i)))/a)) # Normal/Tangential velocity
        v_r.append(vr)
        v_n.append(vn)

    return v_r, v_n

def elements_to_perifocal(ta, a, e, mu, h):
    # Calc perifocal distance in m
   
    r = [] 
    p = []
    q = []
    w = []
    dp = []
    dq = []
    dw = []
   
    
    for theta in ta:
        position = h**2/mu*1/(1+e*np.cos(theta))
        # position = h**2/(mu*(1+(e*np.cos(theta))))

        #perifocal distances
        pn = position*np.cos(theta)
        qn = position*np.sin(theta)
        
        p.append(pn)
        q.append(qn)
        w.append(0)
        # Compute perifocal velocities
        dpn = -np.sin(theta)*mu/h
        dqn = (e+np.cos(theta))*mu/h
        dp.append(dpn)
        dq.append(dqn)
        dw.append(0) 
    
    return p, q, w, dp, dq, dw


def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):
    
    x = []
    y = []
    z = []
    dx = []
    dy = []
    dz = []


    eci = []
    # Compute coordinates
    for k in range(len(p)):

        eci_position = np.array(np.matmul(perifocal_to_eci_matrix(i,raan,argp),[p[k],q[k],w[k]]))
        eci.append(eci_position)
        
        xk = eci_position[0]
        x.append(xk)
        yk = eci_position[1]
        y.append(yk)
        zk = eci_position[2]
        z.append(zk)
    
    for i in dx:
        eci_velocity = list(np.array(np.matmul(perifocal_to_eci_matrix(i,raan,argp),[dp[k],dq[k],dw[k]]))) 
        dx.append(eci_velocity)  

    for i in dy:
        eci_velocity = list(np.array(np.matmul(perifocal_to_eci_matrix(i,raan,argp),[dp[k],dq[k],dw[k]]))) 
        i = i + eci_velocity[1] 
        dy.append(eci_velocity)  
    for i in dz:
        eci_velocity = list(np.array(np.matmul(perifocal_to_eci_matrix(i,raan,argp),[dp[k],dq[k],dw[k]]))) 
        dz.append(eci_velocity)  
    
    
    return x, y, z, dx, dy, dz

def perifocal_to_eci_matrix(i, raan, argp):

    
    # Calculate transformation matrix from perifocal to ECI frame
    r1 = -np.sin(raan)*np.cos(i)*np.sin(argp)+np.cos(raan)*np.cos(argp)
    r2 = -np.sin(raan)*np.cos(i)*np.cos(argp)-np.cos(raan)*np.sin(argp) #note this line had an error, it had np.cos(raan) instead of np.sin(raan)
    r3 = np.sin(raan)*np.sin(i)
    
    r4 = np.cos(raan)*np.cos(i)*np.sin(argp)+np.sin(raan)*np.cos(argp) #basically I think the error is to do with me using the transformation matrix from eci to perifocal instead of the other way around, I got some of the inputs into the matrix mixed up   
    r5 = np.cos(raan)*np.cos(i)*np.cos(argp)-np.sin(raan)*np.sin(argp) #and my brain got mixed up on how matrices (arrays) work in python
    r6 = -np.sin(i)*np.cos(argp)

    r7 = np.sin(i)*np.sin(argp)
    r8 = np.sin(i)*np.cos(argp)
    r9 = np.cos(i)

    p_to_e = np.array([[r1,r2,r3],[r4,r5,r6],[r7,r8,r9]]) 
    
    return p_to_e
