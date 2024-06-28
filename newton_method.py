import numpy as np

def calculate_anomalies(t_sec, ma, e, n, T):
    mean_anomaly = []
    theta = []
    for i in t_sec:
        M_e = 2*np.pi*i/T
        mean_anomaly.append(M_e)
        E_0 = M_e


        # fE = E_0-(e*np.sin(E_0))-M_e
        # f_prime = 1-e*np.cos(E_0)
        # E_n = E-(fE/f_prime)


        n = 100
        # print(E_n-E)
        tol = 1e-8
        for i in range(1,n):
            fE = E_0-(e*np.sin(E_0))-M_e
            f_prime = 1-e*np.cos(E_0)
            E_1 = E_0-(fE/f_prime)
            # E_0 = E_1
            if abs(E_1-E_0) < tol:
                E = E_1
                i = n
            else:
                E_0 = E_1
    

        A = np.tan(E/2)
        B = np.sqrt((1-e)/(1+e))
        ta = 2*np.arctan(A/B)
        theta.append(ta)
    return mean_anomaly, theta
   