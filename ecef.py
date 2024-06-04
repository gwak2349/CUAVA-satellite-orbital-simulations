import numpy as np

def eci_to_ecef(t_0,t_sec,x,y,z):
    rotation_speed = (np.pi*2*(1+1/365.26)/(24*3600))

    # print(rotation_speed)
    angle_list = []
    # position_ecef = []

    x_ecef = []
    y_ecef = []
    z_ecef = []
    r = []
    alpha = []
    delta = []
    for t in range(len(t_sec)):
        angle = rotation_speed*(t_sec[t]-t_0)
        angle_list.append(angle)
        r1 = [np.cos(angle),np.sin(angle),0]
        r2 = [-np.sin(angle),np.cos(angle),0]
        r3 = [0,0,1]
        eci_to_ecef = np.array([r1,r2,r3])
        vector = np.array([x[t],y[t],z[t]])
        p_ecef = np.matmul(eci_to_ecef,vector)
        # p_ecef = [r1[0]*x[t]+r1[1]*y[t]+r1[2]*z[t],
        #           r2[0]*x[t]+r2[1]*y[t]+r2[2]*z[t],
        #           r3[0]*x[t]+r3[1]*y[t]+r3[2]*z[t]]
        # print(p_ecef)
        x2 = p_ecef[0]
        x_ecef.append(x2)
        y2 = p_ecef[1]
        y_ecef.append(y2)
        z2 = p_ecef[2]
        z_ecef.append(z2)
        magnitude = np.sqrt(x2**2+y2**2+z2**2)
        r.append(magnitude)
        l = x2/magnitude #direction cosine in x-direction
        m = y2/magnitude #direction cosine in y-direction
        n = z2/magnitude #direction cosine in z-direction
        delta_i = np.arcsin(n)
        
        delta.append(delta_i)
        
        if m > 0:
            alpha_i = np.arccos(l/np.cos(delta_i))
            alpha.append(alpha_i)

        else:
            alpha_i = -np.arccos(l/np.cos(delta_i))
            alpha.append(alpha_i)

        # position_ecef.append(p_ecef)
    # for k in range(len(angle_list)):
    #     angle_list[k] = angle_list[k]*180/np.pi


    # for k in range(len(alpha)):
    #     alpha[k] = alpha[k]*180/np.pi
    #     delta[k] = delta[k]*180/np.pi
    
    return x_ecef, y_ecef, z_ecef, r, alpha, delta 

# def convert_radec_to_longlat(alpha, delta):

#     sidereal_time = alpha


# #hour + minute/60 + second/3600 = decimal value
#     longitude = 1
#     latitude = delta

#     return longitude, latitude
