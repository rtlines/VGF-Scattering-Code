#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 14:25:34 2022

@author: rtlines
"""
# import statements go here ###################################################
import random 
import numpy as np
import matplotlib.pyplot as plt
import json
#
#
# Function definitions section ################################################

###############################################################################
#
###############################################################################
def read_input_file(workfile):
   # comments is the coment string from vgfin
   # NUSE is the number of dipoles to use
   # wave is the wavelength in relative uints (as long as all lenght units 
   #      are the same it doesn't matter, but use microns)
   # alpha, beta, and gamma are the euler rotation angles
   # psi is the polarization angle (I think)
   # RAD is the particle semimajor axis (use microns)
   # ER and EI are the complext components of the particle permitivity
   # TD is the dipole size
   # R is a NUSE by 3 array of cell locations
   # D is a NUSE array of cell weights
   #
   # Start reading in data
   #print("start read intput file")
   with open(workfile) as f:
        print("read comment")
        read_data = f.readline()            # read a whole line
        com=read_data
        print(com)
        read_data = f.readline()            # read a whole line
        print("NUSE",read_data)
        NUSE = int(read_data)
        #print("NUSE", NUSE)
        read_data = f.readline()            # read a whole line
        split_string = read_data.split(" ")  #parse the line splitting where
        #print(split_string)
        wave = float(split_string[3])
        alpha = float(split_string[10])
        beta = float(split_string[17])
        gamma = float(split_string[24])
        gamma = float(split_string[31])
        RAD = float(split_string[38])
        #print("W,alpha, beta, gamma, psi, RAD", wave,alpha,beta,gamma,psi,RAD)
        read_data = f.readline()            # read a whole line
        split_string = read_data.split(" ")  #parse the line splitting where
        #print("  ")
        #print(split_string)  
        ER = float(split_string[3])
        EI = float(split_string[9])
        TD = float(split_string[15])
        #print("ER, EI, TD", ER, EI, TD)  
        # make a place for the locations of the dipole cells, we need x,y,z components
        # and make a place for the dipole weights
        R = np.zeros((NUSE,3))  
        D = np.zeros(NUSE)
        for i in range(0,NUSE,1):             # loop to input all the theta and phi values
           #print("I = ",i," NUSE = ", NUSE)
           read_data = f.readline()       # read in a whole line
           #print("read_data in loop",read_data)
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
           #print("loop end")
           #still need 
        #print("after loop")        
        read_data = f.readline()            # read a whole line
        #print("after loop",read_data)
        temp = float(read_data)
        #print(temp)
        read_data = f.readline()            # read a whole line
        #print(read_data)
        factor = float(read_data) 
        #print(factor)
   #print("return to main program")
   return NUSE,com,wave,alpha,beta,gamma,psi,RAD,ER,EI,TD,R,D
###############################################################################
#
###############################################################################
def kvector_components(theta, phi):
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
#
#  Test code
#
######
#NK=2                                  # Asking for 2 vectors
#khatN = np.zeros((NK,3))              # array of components of the vectors
#kth = np.zeros(NK)                    # array of theta values
#kph = np.zeros(NK)                    # array of phi values
#kth[0] = np.pi/2                      # test theta values
#kth[1] = np.pi/3
#kph[0] = 0.0                          # test phi values
#kph[1] = np.pi
#
#N = 1                                 # choose the vector to find components
#print(khatN)                          # print to make sure all zeros
#khatN[N] =  kvector_components(kth[N],kph[N])   # find the components
#print(khatN)                          # print to check it worked
#####                                      
# for N = 0 you should get
#    [[1.000000e+00 0.000000e+00 6.123234e-17]
#    [0.000000e+00 0.000000e+00 0.000000e+00]]
#
# for N = 1 you should get 
#     [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00]
#     [-8.66025404e-01  1.06057524e-16  5.00000000e-01]]
#
# End Test code
###############################################################################
#
###############################################################################
# function to find the components of the k-vectors as a group
def kvectors3(kth,kph,NK):                      
    khatN=np.zeros((NK,3))
    for N in range(NK):          
           khatN[N][0]=np.sin(kth[N])*np.cos(kph[N]) 
           khatN[N][1]=np.sin(kth[N])*np.sin(kph[N])
           khatN[N][2]=np.cos(kth[N])
    return khatN
###############################################################################
#
###############################################################################

def move1kv(mcount,khatN,kth,kph,NK,iseed,divd):
    # mcount is the index of a particular k-vector number
    # khatN is the Cartisian components of the k-vector 
    # kth and kph are our angle arrays we input above
    # NK is the number of k-vectors
    # iseed is a randum number seed (do this later)
    # divid is the amount we are dividing pi by to get our increment
    #
    # start by defining an angular increment abount int the theta direction
    #
    THinc = np.pi/divd                 #Theta increment
    PHinc = 2*THinc                    #Phi direction increment is twice as big
    temp1 = THinc*random.uniform(-1,1) # this is our theta delta
    temp2 = PHinc*random.uniform(-1,1) # this is our phi delta
    #
    kth[mcount] = kth[mcount]+ temp1   #Add the delta theta to the partucuar
                                       #  theta value
    kph[mcount] = kph[mcount] + temp2  #Add the delta phi to the partucuar
                                       #   phi value
    # we have a new k-vector angle, but the change could make us get strange
    # vectors. We want positive angles, but we also want to adjust the angles
    # either in the positive or the negative direction. This sometimes could 
    # make negagive vector angles. Check for this and correct it if it happens
   
    if (kth[mcount] > np.pi):
        kth[mcount] = 2*np.pi-kth[mcount]
    if (kph[mcount] > np.pi*2):
        kph[mcount] = kph[mcount]-np.pi*2
    if (kth[mcount] < 0.0):
        kth[mcount] = -kth[mcount]
    if (kph[mcount] < 0.0):
        kph[mcount] = -kph[mcount]
    # Now we need to find the components. A function called kvectors4 does that
    # in the original code. kvectors4(mcount, khatN, kth, kph, NK)
    #
    khatN[mcount]=kvector_components(kth[mcount], kph[mcount])
    #
    #or we could just do it here
    # we should worry about whether we are starting arrays at 0 or 1!!!!! @@@@
    #khatN[mcount,0] = np.sin(kth[mcount])*np.cos(kph[mcount]) 
    #khatN[mcount,1] = np.sin(kth[mcount])*np.sin(kph[mcount])
    #khatN[mcount,2] = np.cos(kth[mcount])
    #
    # and we are done. We have a set of k-vectors with one vector randomly
    # moved.
    return kth, kph, khatN 
###############################################################################
#
###############################################################################

def read_kv(workfile):
   #workfile ="Test_2_5_5.kv"
   with open(workfile) as f:
      read_data = f.readline()            # read a whole line
      #print(read_data)                    # The first line is just the number of kvectors
      NK=int(read_data)                   # turn the whole line into a number
      #print(NK)
#      kth = np.zeros(NK)                  # make arrays for the kvector components
#      kph = np.zeros(NK)
#      khatN = np.zeros((NK,3))              # actual components of the vectors
      for i in range(0,NK,1):             # loop to input all the theta and phi values
          read_data = f.readline()       # read in a whole line
          #print(read_data)
          split_string = read_data.split(" ")  #parse the line splitting where
                                               #there are spaces
                                               #Each line starts with a sppace
                                               #then a number, then a space
                                               #then a number. So we will get
                                               #three parts to our split
                                               #a space, a string that is a 
                                               #number, then a second string 
                                               #that is a number.
          #print(split_string)
          theta=split_string[1]        # ignore the first split piece and 
                                       # copy the first number string into 
                                       # a variable
          phi=split_string[2]          # and do the same for the second number
                                       # string
          #print(theta, phi)     
          kth[i] = float(theta)        # now convert the string into a number
                                       # and put it into the array
          kph[i] = float(phi)          # and do the same for the second string
          #print(kth[i],kph[i])
          
      read_data = f.readline()       # The last line has dummy strings
                                     # for kth and kph but has at the end 
      #print(read_data)
      split_string = read_data.split(" ")
      #print(split_string)
      ERR=float(split_string[3])
      ERRlast=float(split_string[6])
      mcount=int(split_string[21])
      kcount=int(split_string[32])
      #print(ERR, ERRlast, mcount, kcount)
      return NK,kth,kph, ERR, ERRlast, mcount, kcount

# now we have the kth and kph values in arrays just like the original
# code wanted.
###############################################################################
#
###############################################################################
def ROT(alpha, beta, gamma):
    RRR = np.array(
      [[ np.cos(alpha)*np.cos(beta)*np.cos(gamma)-np.sin(alpha)*np.sin(gamma),
         np.sin(alpha)*np.cos(beta)*np.cos(gamma)+np.cos(alpha)*np.sin(gamma),
        -np.sin(beta)*np.cos(gamma)],
       [-np.cos(alpha)*np.cos(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma),
        -np.sin(alpha)*np.cos(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma),
         np.sin(beta)*np.sin(gamma)],
       [ np.cos(alpha)*np.sin(beta),
         np.sin(alpha)*np.sin(beta),
         np.cos(beta)]])
    return RRR
###############################################################################
#
###############################################################################
# Calculates a Dirac Delta Function
def delta(alpha, beta):
    d=0.0
    if ( alpha == beta ):
        d=1.0
    return d  
###############################################################################
#
###############################################################################
# Calculates a double Dirac Delta Function
def dd(alpha, i, beta, j):
    d=0.0
    if ((alpha == beta) and (i == j)):
        d=1.0
    return d

###############################################################################
#
###############################################################################
#    Function to calculate the self term contribution termed GAMMA    
#      in the IBM write-up. The form for GAM  is taken from the work  
#      of B. T. Draine and J. Goodman, Astrophysical Journal, 405:    
#      685-697, 1993 March 10. In my experience, the Draine and   
#      Goodman formulation is the better of the three.  But it is more 
#      complicated.  It looks like this uses Goedecke and O'Brian.
#      Note the sign change due the time dependance in VGF being exp(-iwt). 
#      This different from Goedecke and O'Brien's choice. 
#
#            GAM=(3./(4.*PI))**(2./3.)*(kd)**2 + CI*kd**3/(2.*PI)          
#
#      Note that in the Goedecke and O'Brien series, Goedecke and     
#      O'Brien expand the exponential in the self term integral and   
#      throw away most of the trems.                               
#           real a                                                        
#           complex temp                                                  
#           a=d*(3./(4.*PI))**(1./3.)                                     
#           temp=CI*k*a                                                   
#           GAM=2.*((1.-CI*k*a)*cexp(temp)-1.)      
#      Single Point Test 2/16/2022    
def GAM(d,k,EPS):
       GAM=0.0+0.0j
       kd=k*d
       b1=(3./(4.*np.pi))**(2./3.)
       GAM=b1*(kd)**2 + 1j*kd**3/(2.*np.pi)          
       return GAM
###############################################################################
#
###############################################################################
# Calculate the Dyadic Green's Function for a dipole 
#
def GG(R,k,d,EPS,a,i,b,j):
    # Inputs:
    #    R matrix of dipole locations (Nuse,3)
    #    k the wave number
    #    d is a specific cell weighting (passed in as a single value)
    #    EPS is the complex permitivity
    #    a, i, b, and j are indexes to get the right component of the right 
    #    dopole location.
    # set up the complex temporary variables and zero them
    PHZ = 0.0+0.0j
    t1 = 0.0+0.0j
    t2 = 0.0+0.0j
    temp = 0.0+0.0j
    ci = 0.0+1.0j

    #
    K2=k*k 
    # make a space for dipole displacement vectors 
    Rab = np.zeros(3)
    Rhat = np.zeros(3)
    if b != a:  
        # @@@ Check for error in this branch
        # calculat4e the separation distance between two dipoles
        #   Rmn = Rn - Rm and RMAG = |Rmn|
        Rab[0] = R[a][0] - R[b][0]
        Rab[1] = R[a][1] - R[b][1]
        Rab[2] = R[a][2] - R[b][2]
        # now find the magnitude
        RMAG = Rab[0]**2 + Rab[1]**2 +Rab[2]**2 
        RMAG = RMAG**0.5
        # make a unit vector in the Rmn direction
        Rhat[0] = Rab[0]/RMAG
        Rhat[1] = Rab[1]/RMAG
        Rhat[2] = Rab[2]/RMAG
        # our exponent is i(k dot r) form the exponent
        temp = ci*k*RMAG
        PHZ = np.exp(temp)
        #
        t1 = (K2/RMAG)*(delta(i,j) - Rhat[i]*Rhat[j])
        t2 = (ci*k/RMAG**2 - 1.0/RMAG**3) * (delta(i,j)-3.0*Rhat[i]*Rhat[j])
        G = PHZ * (t1 + t2)
    else:    
        G = 4.0*np.pi*GAM(d,k,EPS)/(3.0*d**3)
    return G
#
###############################################################################
#
###############################################################################
# Function to calculate the trial funciton expansion functions 
#    CPSI=cexp(i*k*khatN.R(b))  
# 
# Single point checked 02/16/2022 

def CPSI(R,khatN,k,mm,N,b) :
    KDR=khatN[N].dot(R[b])
    temp = 1j*mm*k*KDR
    return np.exp(temp)
    
     

###############################################################################
#
###############################################################################
# function to print T1 to a file
def print_T1_to_file(T1,NUSE, NK):
    print("Printing T1 data to T1_file.txt")
    with open('T1_file.txt','w') as Tf:
       for a in range(NUSE):
            for i in range(3):
                for N in range (NK):
                    for j in range(3):
                        temp = str(a)+", "+str(i)+", "+str(N)+", "+str(j)
                        temp = temp + ", " + str(T1[a][i][N][j])+"\n" 
                        Tf.write(temp)            # write a whole line
                        # End of the j loop         
                    # End of the N loop
            # End of the a loopallowed
       # End of file open
###############################################################################
#
###############################################################################
# function to print Y to a file
def print_Y_to_file(Y, NK):    
    print('Printing Y data to Y_file.txt')
    with open('Y_file.txt','w') as Yf:
        for M in range(NK):
              for l in range(3):  
                    temp = str(M) +", "+str(l)+", "
                    temp = temp + str(Y[M][l])+"\n"
                    Yf.write(temp)            # write a whole line
                    # End l loop
              # End M loop   
        # Close file
###############################################################################
#
###############################################################################
# function to print H to a file    
def print_H_to_file(H,NK):      
    print('Printing H data to H_file.txt')
    with open('H_file.txt','w') as Hf:
         for M in range (NK):
             for l in range (3):
                 for N in range(NK):
                     for i in range(3):
                         temp = str(M)+", "+str(l)+", "+str(N)+", "+str(i)
                         temp = temp + str(H[M][l][N][i])+"\n"
                         Hf.write(temp)            # write a whole line
                         # End i loop
                     # End N loop
                 # End l looop
             # End M loop
         # Close File    
###############################################################################
#
###############################################################################
# function to print T1 to a file
def print_GG_to_file(R, k, eps):
    print("Printing GG data to GG_file.txt")
    dtemp = 1.24
    with open('GG_file.txt','w') as Gf:
       for a in range(NUSE):
            for i in range(3):
                for b in range (NK):
                    for j in range(3):
                        temp = str(a)+", "+str(i)+", "+str(b)+", "+str(j)
                        temp = temp + ", " + str(GG(R,k,dtemp,eps,a,i,b,j))+"\n" 
                        Gf.write(temp)            # write a whole line
                        # End of the j loop         
                    # End of the N loop
            # End of the a loop
       # End of file open
###############################################################################
#
###############################################################################
# function to print bb to a file
def print_bb_to_file(bb, NK):
    print("Printing bb data to bb_file.txt")
    with open('bb_file.txt','w') as bf:
       for i in range(NK*3):
            temp = str(i)+", "+str(bb[i])+"\n" 
            bf.write(temp)            # write a whole line
            # End of the i loop
       # End of file open
###############################################################################
#
###############################################################################
# function to print aa to a file
def print_aa_to_file(aa, NK):
    print("Printing aa data to aa_file.txt")
    with open('aa_file.txt','w') as aaf:
       for i in range(NK*3):
           for j in range(NK*3):
               temp = str(i)+", "+str(j)+", "+ str(aa[i][j])+"\n" 
               aaf.write(temp)            # write a whole line
            # End of the i loop
       # End of file open       
###############################################################################
#
###############################################################################
# function to tell if our field expansion coefficients are good
#   so we can end the Monty Carlo loop
def Calculate_PHI(An, H, Y, E0, NK, NUSE):
    CPHI=0.0+0.0j
    for N in range(NK):
        for j in range(3):
            for M in range(NK):
                for l in range (3):
                    CPHI = CPHI + np.conjugate(An[M][l]*H[M][l][N][j]*An[N][j])
                    # End l loop
                # End M loop
            # End j loop
        # End N loop
    for j in range (3):
        for N in range (NK):
            CPHI = CPHI - (np.conjugate(Y[N][j])*An[N][j]                     \
                           +Y[N][j]*np.conjugate(An[N][j]))
            # End N loop
        # End j loop
    for a in range(NUSE):
        for i in range(3):
            CPHI - CPHI + E0[a][i]*np.conjugate(E0)[a][i]
            # End i loop
        # End a loop        
    return CPHI
    #
###############################################################################
#
###############################################################################

def ecalc (NUSE,R,E, An,khatN,X,NK,K,mm)  :
    # NUSE is the number of dipoles we used to repdresent the particle
    # R is the array of dipole locations
    # E is the array of electric field vectors
    # F is the array of F-field vectors
    # An is the array of amplitudes for our field expansion
    # khatN is the array of k vector components
    # X is the electric suseptibility
    # NK is the number of k vectors
    # mm I think is the complex index of refraction
    #
    # Form the F-field vectors
    F = np.zeros((NUSE,3),dtype = complex)
    for b in range(0,NUSE):
        for j in range(3):
            F[b][j]=(0.0,0.0)
            for N in range(0,NK):
                F[b][j]=F[b][j]+An[N][j]*CPSI(R,khatN,k,mm,N,b)          
    #
    #Done, Calculate the internal E-field to output it.
    #
    for m in range(0, NUSE):
        for I in range(3):
            E[m][I]=F[m][I]/(1.+(4.*np.pi)*X/3.)
    #
    return E
###############################################################################
#
###############################################################################
# Caclulate the T1 matrix
#  Checked 02/18/2022
def Calculate_T1_Matrix( D, R, k, khatN, eps, mm, W, NUSE, NK ):
    for a in range(NUSE):
        for i in range(3):
            for N in range (NK):
                  for j in range(3):
                      T1[a][i][N][j] = 0.0+0.0j
                      for b in range (NUSE):
                          dtemp = D[b]
                          T1[a][i][N][j] =  T1[a][i][N][j] +                  \
                                        ( dd(a,i,b,j)-D[b]**3 * W             \
                                        * GG(R,k, dtemp,eps,a,i,b,j))         \
                                        * CPSI(R,khatN, k, mm, N, b)
                          # End of the b loop
                          
                          #print ("T1 loop ",a,i,N,j,T1[a][i][N][j])
                      # End of the j loop         
                  # End of the N loop
            # End of the i loop
        # End of the a loop
    return T1
###############################################################################
#
###############################################################################



# Main function goes here
#
# Test reading in the k-vectors
#
# array of theta and pi values to read in.  These are the direction angles
#   for the unit k-vectors.
#
# Define Variables
#
# Array Limits:  First let's set limits on how big the arrays can be
KMAX = 2000     # maximum number of k-vectors
NMAX =2000      # maximum number of dipoles
MMAX = 2        # maximum number of Monte Carlo itterations (starts at 1?)
#
# Convinient consDEG=PI/180.0tants
DEG=np.pi/180.0     # conversion from degress to radians
#
####### read in the dipole locatons and data from the input file
#
# The code needsto input some things.  I will hard code them here for now @@@@@
ERR0 = 0.05
divd = 100
#
# Now read the input file.  The file name should be an input, but hard code
#    for now @@@@@
#
print("read input file")
comments = "Not Set"          #is the coment string from vgfin
NUSE = 0.0              # the number of dipoles to use
wave = 0.0              # the wavelength in relative uints (as long as all lenght units 
                   #    are the same it doesn't matter, but use microns)
alpha = 0.0            # The tree Euler angles describing the particle orientation
beta = 0.0
gamma = 0.0 
psi = 0.0                # the polarization angle (I think)
RAD = 0.0                # is the particle semimajor axis (use microns)
ER = 0.0                 # the complext components of the particle permitivity
EI = 0.0
TD = 0.0                 # the dipole size
R = np.zeros((NMAX,3))    # a NUSE bcom,wave,alpha,beta,gamma,psi,RAD,ER,EI,TD,R,Dy 3 array of cell locations
D = np.zeros(NMAX)         # a NUSE array of cell weights
workfile = "Figure3.in" # @@@ needs to be an input
[NUSE,comments,wave,alpha,beta,gamma,psi,RAD,ER,EI,TD,R,D] = read_input_file(workfile)
# Test to see if it worked
#print(comments)
#print(NUSE)
#print(wave, alpha, beta, gamma, psi, RAD)
#print(ER, EI, TD)
#print(R)
#print(D)
#
#
###### Read in the K-Vectors
print('read in the k-vector file')
NK = 0                                # number of k-vectors
ERR = 0                               # Monte Carlo Loop Test Value (MCTV)
ERRlast = ERR                         # Last value of MCTV to see if we are 
                                      #   converging to a better answer
kth = np.zeros(KMAX)                  # make arrays for the kvector components
kph = np.zeros(KMAX)
# array of k-vector cartisian components.
khatN = np.zeros((KMAX,3))              # arrat of components of the vectors
kworkfile ="Figure3.kv"   #@@@ needs to be an input
[NK,kth,kph, ERR, ERRlast, mcount, kcount] = read_kv(kworkfile)       # Function call to get the k-vectors
# Test to see if it worked
#print(NK)
#for i in range(NK):
#    print(kth[i],kph[i])
#
# find the components of the k-vectors
khatN = kvectors3(kth,kph,NK) 
###### Now we have the data input, start calculations
#
## Find the index of refraction, epsilon, and the suseptibility, X
eps = complex(ER,EI)
#print(eps)
mm = np.sqrt(eps)
#print(mm)
X = complex(0.0,0.0)
X = (eps - 1)/(4.0*np.pi)
#print(X)
#
###### Begin the incident field calculations
#
k=2.0*np.pi/wave         # wave number
#
print('rotate the input field into particle frame, calc field at each dipole')
# covert angles from degrees to radians
alpah = alpha * DEG
beta = beta * DEG
gamma = gamma * DEG
psi = psi * DEG
#
# set up the rotation matrix with a call to ROT
#
RRR = ROT(alpha, beta, gamma)
#print(RRR)
#
# K-hat is in the z direction in the lab frame, we need to rotate it into the
#   particle frame
V=np.zeros(3)                                 # temporary storage
V[0] = 0.0                                    
V[1] = 0.0
V[2] = 1.0
khat = RRR.dot(V)
# Thus E0hat must be in the x-y plane in the lab. Rotate it into the particle
#   frame     
V[0] = np.cos(psi)
V[1] = np.sin(psi)
V[2] = 0.0
E0hat = RRR.dot(V)
print("E0hat ",psi, E0hat)
#
# Calculate the W factor    Wcalc=X/(1.+(4.*PI/3.)*X)
W = X/(1.0+(4.*np.pi/3.)*X)
#
###### Calcuate inital and incident field
#
# We need a place to put the calcualted electric field at each dipole location
#
E0 = np.zeros((NUSE,3), dtype=complex)
#
# we need some complex temporary varuables
temp = 0+0j
C = 0+0j
#
# first calcualte the r dot khat part of the exponential
#
for i in range (NUSE):
    RDK = np.dot(khat, R[i])       # khat dot R
    #print("RDK ", RDK)
    temp = 1j*k*RDK                # i times the wave number times RDK
    #print("temp ",temp)
    C = np.exp(temp)               # take the exponential
    #print ("C ",C)
    E0[i] = C*E0hat                # it is a plane wave, so we need the 
                                   # complex amplitude multiplied by the 
                                   # plane wave exponential
    #print('E0',E0[i],'\n')
#
# We are going to calcualte the field in a big loop
#    Set up the big loop.
#    We will need spaces for the large matracies to go into
#
T1 = np.zeros((NUSE,3,NK,3),dtype = complex)
Y  = np.zeros((NK,3), dtype = complex)
H  = np.zeros((NK,3,NK,3), dtype = complex)
An = np.zeros((NK,3), dtype = complex)
bb = np.zeros((3*NK), dtype = complex)
aa = np.zeros((3*NK, 3*NK), dtype = complex)
xx = np.zeros((3*NK), dtype = complex)
#
# Code to check the functions we will use, comment out later
# Test CPSI ################################################
#print('testing CPSI')
#N=1
#b=1
#testCPSI = CPSI(R,khatN,k,mm,N,b)
#print('khatN, Rb, mm, k', khatN[N],R[b],mm,k )
#print('CPSI test = ',testCPSI)
#print('Calculationg internal fields...')
# Test dd ################################################
#print('Testing dd')
#a=0
#i=1
#b=0
#j=1
#print(a, i, b, j, 'dd= ', dd(a,i,b,j))
# Test GAM ################################################ 
#print('Testing GAM')
#dtemp = 1.24
#print('d, k, eps ',dtemp,k,eps)
#print('GAM = ',GAM(dtemp,k,eps))

# Test GG ################################################   
#print('Testing GG')
#a=0
#i=1
#b=2
#j=0
#dtemp = 1.24
#tempGG = GG(R,k,dtemp,eps,a,i,b,j)
#print(a,i,b,j,dtemp,' GG= ',tempGG)
#print_GG_to_file(R, k, eps)

while (ERR > ERR0) and (mcount < MMAX):

    print('Monte Carlo Loop mcount = ', mcount)
    #ERR came from the k-vector file. We read it in before
    #kcount is the number of kvector tries to loop over. ikcount seems to have 
    #   allowed us to start not at the first k-vector.  I don't see why we would
    #   do that so let's try without it. The kcount loop gives us the number
    #   of random tries for best fit for the field.  But if we get a good one
    #   we drop out of the loop.
    for kcount in range (1):        # should be NK @@@@ put it back after testing T1
        ##################### Calculate the T1 matrix -------------------------
        print('Calculating T1, NK =',NK,' kcount =', kcount)
        T1 = Calculate_T1_Matrix( D, R, k, khatN, eps, mm, W, NUSE, NK )
        # output the T1 matrix for debugging
        print_T1_to_file(T1,NUSE, NK)
#        for a in range(NUSE):
#            for i in range(3):
#                for N in range (NK):
#                    for j in range(3):
#                        T1[a][i][N][j] = 0.0+0.0j
#                        for b in range (NUSE):
#                            dtemp = D[b]
#                            T1[a][i][N][j] =  T1[a][i][N][j] +                \
#                                ( dd(a,i,b,j)-D[b]**3 * W                     \
#                                 * GG(R,k, dtemp,eps,a,i,b,j))                \
#                                 * CPSI(R,khatN, k, mm, N, b)
                            # End of the b loop
                        #print ("T1 loop ",a,i,N,j,T1[a][i][N][j])
                        # End of the j loop         
                    # End of the N loop
        # End of the a loop
        ################# Calculate the Y matrix ------------------------------
        # Calculate the H and Y matricies
        print('Building the Y matricx')
        for M in range(NK):
              for l in range(3):
                    Y[M][l] = 0.0+0.0j    
                    for a in range(NUSE):
                        for i in range(3):
                            Y[M][l] = Y[M][l] + np.conjugate(T1[a][i][M][l])  \
                                      *E0[a][i]
                            #print("M, l, a, i",M, l,a,i, T1[a][i][M][l],E0[a][i]  )
                            # End i loop
                        # End a loop
                    #print("M, l, Y[M][l]",M, l, T1[a][i][M][l],E0[a][i]  )
                    # End l loop
                # End M loop   
        print_Y_to_file(Y, NK)
        # Y matrix checked on 02/18/2022
        ################# Calculate the Y matrix ------------------------------
        print('Building the H matrix')
        for M in range (NK):
              for l in range (3):
                    for N in range(NK):
                        for j in range(3):
                            H[M][l][N][j] = 0.0+0.0j
                            for a in range(NUSE):
                                for i in range(3):
                                    H[M][l][N][j] =  H[M][l][N][j]            \
                                        + np.conjugate(T1[a][i][M][l])        \
                                        * (T1[a][i][N][j])
                                    # End i loop
                                # End a loop
                            # End j loop
                        # End N loop
                    # End l looop
              # End M loop
        print_H_to_file(H,NK)
        ################# Now get ready for the matrix inversion---------------
        # Now get ready for the matrix inversion
        # But we need H in a good form for the matrix inverter
        # Try to reform it into a NK*3 by NK*3 matrix
        print('reforming the matrix for solving ')
        for n in range(NK):
              for i in range(3):
                  m = 3*(n-0)+i     # in fortran arrays start with 1, but 
                                    # our python arrays start with 0
                  bb[m] = Y[n][i]
                  #print('m, i, m, bb ',n, i, m, bb[m])
                  for nnp in range (NK):
                      for j in range (3):
                          mp = 3*(N-0)+j
                          aa[m][mp] = H[n][i][nnp][j]
                          #print('m, mp, aa[m][mp]',m, mp, aa[m][mp])
                          # End j loop
                      # End np loop
                  # End i loop
              # End n loop
        print_bb_to_file(bb, NK)
        print_aa_to_file(aa, NK)
        # Now we want to solve the matrix equation aa * xx = bb
        # python is supposed to be able to do this with its linear algebra
        # function linalg.solve(aa,bb)
 #       print('matrix solve')
#        xx = np.linalg.solve(aa,bb)
#        # now take our solutino and put it back into matrix component format
#        for N in range (NK):
#                for i in range (3):
#                    M = 3*(N-0)+i
#                    An[n][i] = xx[M]
#                    # End i loop
#                # End N loop
#            # Now we need to test our solution to see how good it is
#            # 
#            PHI = CPHI(An, H, Y, E0, NK, NUSE)  
#            ERR = PHI * np.conjugate(PHI)
#            if (ERR>ERR0):      # Not done with the MC loop, keep going
#               if (ERR > ERRlast):  # Reject: our try isn't so good, reset the k-vectors
#                  ERRlast = ERR     # update the ERR history
#                  [NK,kth,kph, ERR, ERRlast, mcount, kcount] = read_kv(kworkfile)       # Function call to get the k-vectors
#                  for N in range(NK):
#                      khatN[N]=kvector_components(kth[N], kph[N])
#               else:               # Acept:  The try was better, keep it
#                  ERRlast = ERR    # update the ERR history
#                  # because this is an accept, calcualte the fields
#                  E = ecalc (NUSE,R,E, An,khatN,X,NK,K,mm) 
#                  #If we want more Monty carlo tries, keep going, change the k-vectors
#                  if (kcount < NK):
#                      [kth, kph, khatN]=move1kv(mcount,khatN,kth,kph,NK,iseed,divd)
#                 
#            else:               # Accept the MC try
#               E = ecalc (NUSE,R,E, An,khatN,X,NK,K,mm)   
#               break # to exit the kcount loop because we feel we are done
#        #
#        print('end of the kcount loop')        
        # End of kcount loop                        
    #
    # End the Monte Carlo loop
    mcount = mcount +1  
    print ('End Monte Carlo Loop, mcount =',mcount, 'ERR =',ERR)              
# And that is the end of the program
# Of couse I haven't saved off the E-Fields yet.  So there is not output file
#  yet.

   
   







