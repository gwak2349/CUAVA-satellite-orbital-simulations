import numpy as np
import scipy as sp


def compute_anomalies(t_sec, ma, e, n):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = n*t_sec
    
    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return E-e*np.sin(E)
    #Defining the derivative of Kepler's equation
    def derivative_kepler_equation(E):
        return 1-e*np.cos(E)
    
    init_guess = ma
    
    fE = kepler_equation(init_guess)-ma #initial function
   
    f_prime = derivative_kepler_equation(init_guess) #derivative of the initial function
    
    ta_list = []
    ea_list = []

    for j in ma_vec:

        E_0 = init_guess

        #Newton's Method
        fE = E_0-(e*np.sin(E_0))-j
        
        f_prime = 1-e*np.cos(E_0)
        m = 100
        # print(E_n-E)
        tol = 1e-8
        for i in range(1,m):
            E_1 = E_0-(fE/f_prime)
            E_0 = E_1
            fE = E_0-(e*np.sin(E_0))-j
            f_prime = 1-e*np.cos(E_0)
            
            if abs(E_1-E_0) < tol:
                E = E_1
                ea_list.append(E)
            else:
                E_1 = E_0
          
        #True anomaly
        # print(E)
        A = np.tan(E/2)
        B = np.sqrt((1-e)/(1+e))
        
        theta = 2*np.arctan(A/B)
        ta_list.append(theta)

      
    

    # Wrap anomalies to -pi:pi
    for i in ta_list:
        i = float(np.arctan2(np.sin(i), np.cos(i)))
        

    for i in ea_list:
        i = float(np.arctan2(np.sin(i), np.cos(i)))
        

    for i in ma_vec:
        i = float(np.arctan2(np.sin(i), np.cos(i)))

    ea = ea_list 
    ta = ta_list
    
        
    

    return ta, ea, ma_vec


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
        dpn = np.sin(theta)*mu/h
        dqn = (e+np.cos(theta))*mu/h
        dp.append(dpn)
        dq.append(dqn)
        dw.append(0) 
    
    

    return p, q, w, dp, dq, dw


def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):
    

    # Transform coordinates to ECI frame
    # x = np.zeros_like(p)
    # y = np.zeros_like(p)
    # z = np.zeros_like(p)

    # dx = np.zeros_like(p)
    # dy = np.zeros_like(p)
    # dz = np.zeros_like(p)
    
    x = []
    y = []
    z = []
    dx = []
    dy = []
    dz = []


    eci = []
    # Compute coordinates
    for k in range(len(p)):
        # print(perifocal_to_eci_matrix(i,raan,argp))
        # print([p[k],q[k],w[k]])
        eci_position = np.array(np.matmul( perifocal_to_eci_matrix(i,raan,argp),[p[k],q[k],w[k]]))
        eci.append(eci_position)
        
        xk = eci_position[0]
        x.append(xk)
        yk = eci_position[1]
        y.append(yk)
        zk = eci_position[2]
        z.append(zk)
    
    for i in dx:
        eci_velocity = list(np.array(np.matmul( perifocal_to_eci_matrix(i,raan,argp),[dp[k],dq[k],dw[k]]))) 
        dx.append(eci_velocity)  

    for i in dy:
        eci_velocity = list(np.array(np.matmul( perifocal_to_eci_matrix(i,raan,argp),[dp[k],dq[k],dw[k]]))) 
        i = i + eci_velocity[1] 
        dy.append(eci_velocity)  
    for i in dz:
        eci_velocity = list(np.array(np.matmul( perifocal_to_eci_matrix(i,raan,argp),[dp[k],dq[k],dw[k]]))) 
        dz.append(eci_velocity)  
    
    
    return x, y, z, dx, dy, dz

def perifocal_to_eci_matrix(i, raan, argp):

    
    # Calculate transformation matrix from perifocal to ECI frame
    r1 = -np.sin(raan)*np.cos(i)*np.sin(argp)+np.cos(raan)*np.cos(argp)
    r2 = -np.sin(raan)*np.cos(i)*np.cos(argp)-np.cos(raan)*np.sin(argp)
    r3 = np.sin(raan)*np.sin(i)

    r4 = np.cos(raan)*np.cos(i)*np.sin(argp)+np.cos(raan)*np.cos(argp)
    r5 = np.cos(raan)*np.cos(i)*np.cos(argp)-np.sin(raan)*np.sin(argp)
    r6 = -np.cos(raan)*np.sin(i)

    r7 = np.sin(i)*np.sin(argp)
    r8 = np.sin(i)*np.cos(argp)
    r9 = np.cos(i)

    
    p_to_e = np.array([[r1,r2,r3],[r4,r5,r6],[r7,r8,r9]]) 
    
    
    return p_to_e
