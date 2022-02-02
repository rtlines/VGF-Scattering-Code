import numpy as np

workfile = "test_output_file.out"

with open(workfile) as f:
     read_data = f.readline()            # read a whole line
     print(read_data)
     read_data = f.readline()            # read a whole line
     print(read_data)
     NUSE = int(read_data)
     print(NUSE)
     read_data = f.readline()            # read a whole line
     split_string = read_data.split(" ")  #parse the line splitting where
     print(split_string)
     W = float(split_string[3])
     alpha = float(split_string[10])
     beta = float(split_string[17])
     gamma = float(split_string[31])
     psi = float(split_string[38])
     print(W,alpha,beta,gamma,psi)
     read_data = f.readline()            # read a whole line
     split_string = read_data.split(" ")  #parse the line splitting where
     print(split_string)  
     ER = float(split_string[3])
     EI = float(split_string[9])
     TD = float(split_string[16])
     print(ER, EI, TD)
     # make a place for the locations of the dipole cells, we need x,y,z components
     # and make a place for the dipole weights
     R = np.zeros((NUSE,3))
     D = np.zeros(NUSE)
     for i in range(0,NUSE,1):             # loop to input all the theta and phi values
          read_data = f.readline()       # read in a whole line
          #print(read_data)
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
          print(split_string)
          R[i][0] = float(split_string[3])
          R[i][1] = float(split_string[5])
          R[i][2] = float(split_string[8])
          D[i] = float(split_string[11])
          #still need 
          
          
