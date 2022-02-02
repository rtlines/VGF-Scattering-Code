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
#
#
# Function definitions section ################################################
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
    khatN[mcount]=kvectors_components(kth[mcount], kph[mcount])
    #or we could just do it here
    # we should worry about whether we are starting arrays at 0 or 1!!!!! @@@@
    #khatN[mcount,0] = np.sin(kth[mcount])*np.cos(kph[mcount]) 
    #khatN[mcount,1] = np.sin(kth[mcount])*np.sin(kph[mcount])
    #khatN[mcount,2] = np.cos(kth[mcount])
    #
    # and we are done. We have a set of k-vectors with one vector randomly
    # moved.
###############################################################################
#
###############################################################################
def read_kv(kth,kph, workfile):
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
      ERR=int(split_string[7])
      ERRlast=int(split_string[14])
      mcount=int(split_string[25])
      kcount=int(split_string[36])
      #print(ERR, ERRlast, mcount, kcount)
      return NK

# now we have the kth and kph values in arrays just like the original
# code wanted.
###############################################################################
#
###############################################################################
#
#
#
# Main function goes here
#
# Test reading in the k-vectors
#
# array of theta and pi values to read in.  These are the direction angles
#   for the unit k-vectors.
KMAX = 2000
NK=0
kth = np.zeros(KMAX)                  # make arrays for the kvector components
kph = np.zeros(KMAX)
# array of k-vector cartisian components.
khatN = np.zeros((KMAX,3))              # actual components of the vectors
workfile ="Test_2_5_5.kv"
NK = read_kv(kth,kph, workfile)
print(NK)
for i in range(NK):
   print(kth[i],kph[i])




