import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

def eci_plot(x,y,z,R_earth):

    fig = plt.figure() 
    ax = fig.add_subplot(projection='3d')

    #setting up coordinates for the Earth
    u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
    x_Earth = R_earth*np.cos(u) * np.sin(v)
    y_Earth = R_earth*np.sin(u) * np.sin(v)
    z_Earth = R_earth*np.cos(v)

    ax.plot(x,y,z, color="green") #plots the orbit in eci coordinates
    ax.plot_surface(x_Earth, y_Earth, z_Earth, color="blue") #plots the spherical Earth
    ax.set_aspect('equal')
    plt.show()


