import math
import re
import read_tle as rt

def sidereal_time(epoch_year):
    filename = "ESA_XMM_tle.txt"
    line_1, line_2 = rt.read_tle(filename)

    t_epoch = line_1[3]


    t_split = re.split('([^a-zA-Z0-9])',t_epoch) #splitting the epoch year from the epoch integer day

    time = t_split[2]
    ty1 = re.split("",time)
    bp1 = ty1[1:9]



    angle_decimal = float(''.join(bp1))*24*1e-8
    epoch_year = re.split("",t_split[0])

    year = epoch_year[1:3]
    day = epoch_year[3:6]
    
    nu = 20
    zu = float(''.join(year))
    
    year1 = float(str(nu)+str(zu))
    m = 0
    d = float(''.join(day))

    UT = angle_decimal

    
    
    J0 = 367*year1-math.floor(7*(year1+math.floor((m+9)/12))/4)+math.floor(275*m/9)+d+1721013.5
    
    # JD = J0+UT/24
    T0 = (J0-2451545)/36525
    theta_G_0 = 100.4606184+36000.77004*T0 + 0.000387933*T0**2 - 2.583*(1e-8)*T0**3
    # print(theta_G_0)
    
    t1 = math.floor(theta_G_0/360)
    theta = theta_G_0-t1*360
    # theta_G = theta_G_0 + 360.98564724*UT/24
    theta_G = theta + 360.98564724*UT/24
    # print(theta_G)
    # print(theta_G_0)
    if theta_G > 360:
        theta_G = theta_G - 360
    # print(J0)
    # print(JD)
    sidereal_time = theta_G*3600
    return sidereal_time