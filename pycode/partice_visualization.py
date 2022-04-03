import numpy as np

def vector_components(theta, phi):
    # Function to find the components of a unit vector given in three dimensions
    #   given the sherical componets angles theta and phi
    #
    # We are going to call this using our kth and kph arrars and put the result
    #    into khatN.  So the call shoudl be something like
    #    khatN(N)=kvectors_components(kth[N], kph[N])
    #
    khat=np.zeros(3)
    khat[0] = np.sin(theta)*np.cos(phi) 
    khat[1] = np.sin(theta)*np.sin(phi)
    khat[2] = np.cos(theta)
    return khat
    # and we are done making the Nth vector components
    
#workfile = "test_output_file.out"
workfile = "Figure3.in"

with open(workfile) as f:
     read_data = f.readline()            # read a whole line
#     print(read_data)
     read_data = f.readline()            # read a whole line
#     print(read_data)
     NUSE = int(read_data)
     #print(NUSE)
     read_data = f.readline()            # read a whole line
     split_string = read_data.split(" ")  #parse the line splitting where
     #print(split_string)
     W = float(split_string[3])
     alpha = float(split_string[10])
     beta = float(split_string[17])
     gamma = float(split_string[31])
     psi = float(split_string[38])
     #print(W,alpha,beta,gamma,psi)
     read_data = f.readline()            # read a whole line
     split_string = read_data.split(" ")  #parse the line splitting where
     #print("  ")
     #print(split_string)  
     #ER = float(split_string[3])
     #EI = float(split_string[9])
     #TD = float(split_string[15])
     #print(ER, EI, TD)
     # make a place for the locations of the dipole cells, we need x,y,z components
     # and make a place for the dipole weights
     R = np.zeros((NUSE,3))
     D = np.zeros(NUSE)
     for i in range(0,NUSE,1):             # loop to input all the theta and phi values
          read_data = f.readline()       # read in a whole line
          #print("read_data",read_data)
          while "  " in read_data:        # make double spaces into sinble spaces
              read_data = read_data.replace("  "," ")
          read_data.lstrip()           # remove leading space
          read_data.lstrip()           # remove leading space
          #print("read_data",read_data)
          split_string = read_data.split(" ")  #parse the line splitting where
                                               #there are spaces
                                               #Each line starts with a space
                                               #then a number, then a space
                                               #then a vextor and then a weight.
                                               # So we will get
                                               # severak parts to our split
                                               #a space, a string that is a 
                                               #number, then a second string 
                                               #that is a number, etc. 
          #print("split_string1",split_string)
          R[i][0] = float(split_string[1])
          R[i][1] = float(split_string[2])
          R[i][2] = float(split_string[3])
          D[i] = float(split_string[4])
          #still need 
     read_data = f.readline()            # read a whole line
     #print(read_data)
     temp = float(read_data)
     #print(temp)
     read_data = f.readline()            # read a whole line
     #print(read_data)
     factor = float(read_data) 
     #print(factor)
     
          
particle_list=[]

from vpython import *
for i in range(NUSE):
    ball = sphere(pos=vector(R[i][0],R[i][1],R[i][2]),radius=D[i]/2)
    particle_list.append(ball)

scene.background=color.white
particle=compound(particle_list)
vect = np.zeros(3)
Deg = np.pi/180.0
theta = 45*Deg    
phi = 90*Deg
vect=vector_components(theta, phi)
print (vect)
#particle.axis = vector(vect[0],vect[1],vect[2])
particle.rotate(angle=93, axis=vec(0,1,0))

