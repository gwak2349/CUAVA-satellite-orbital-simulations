def read_tle(filename):
    txt = open(filename, "r") #opening the TLE text file
    contents = txt.read() #reading the TLE text file
    TLE = []

    for x in contents.split("\n"):
        TLE.append(x.split(" "[-1])) #splitting each element in the TLE  
        
    line_0 = TLE[0]
    line1 = TLE[1]
    line2 = TLE[2]

    line_1 = [x for x in line1 if x != '']
    line_2 = [x for x in line2 if x != '']
    
    return line_1, line_2