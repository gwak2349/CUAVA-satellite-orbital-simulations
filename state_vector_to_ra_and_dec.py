import numpy as np

def state_vector_to_ra_and_dec(x,y,z):
    
    magnitude = []
    delta = []
    alpha = []
    
    
    
    for i in range(len(x)):
        r = np.sqrt(x**2+y**2+z**2)
        r.append(magnitude)
        l = x/magnitude #direction cosine in x-direction
        m = y/magnitude #direction cosine in y-direction
        n = z/magnitude #direction cosine in z-direction
        delta_i = np.arcsin(n)
        
        delta.append(delta_i)
        
        if m > 0:
            alpha_i = np.arccos(l/np.cos(delta_i))
            alpha.append(alpha_i)

        else:
            alpha_i = -np.arccos(l/np.cos(delta_i))
            alpha.append(alpha_i)
            
            
    return alpha, delta

        # position_ecef.append(p_ecef)
    # for k in range(len(angle_list)):
    #     angle_list[k] = angle_list[k]*180/np.pi
        