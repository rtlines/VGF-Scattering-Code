# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 15:25:27 2021

@author: rtlines
"""
import json
import numpy as np

NSID=10

# Devine some arrays to try to stack together to make a 2D array
A0=np.array([[0,0,0]])  
A1=np.array([[1, 2, 3]])
A2=np.array([[7, 8, 9]])

#Now define an empty array called "P"
P=np.empty((0,3),float)
print("P",P)                  # See if P esists by printing it
P=np.append(P,A0, axis=0)     # Append in A0 and see if it shows up in P
print("P",P)                  # Print P to see if it now has A0
P=np.append(P,A1, axis=0)     # Append in A1 and see if it shows up in P
print("P",P)                  # Print P to see if it now has A0 and A1
P=np.append(P,A2, axis=0)     # Append in A2 and see if it shows up in P
print("P",P)                  # Print P to see if it now has A0, A1, and A2

# If that worked, try it in a loop
print("loop")
R=np.empty((0,3),float)       # Define an empty array called "R"
i=1                           # Define the loop index
while i<NSID:                 # Start the loop
     if i<5:                  # Use some criteria to tell if we include data in our array
       temp_array=np.array([[i,i,i]])      # here is some dummy data to putin
                                           # for us it will be our x, y and z values
                                           # instead of i,i,i we would have
                                           # x,y,z
       R=np.append(R,temp_array, axis=0)   # Now append our temporary array
                                           # into R
     i=i+1                                 # increment the loop counter


print("R",R)


# I don't think we need the rest of this, I just wanted stuff to put into 
# the jason file
AnArray=np.array([10, 11,12,13])
AList1=[1,2,3]
AList2=[4,5,6]

vgfin_data=[{'wavelength':390, 
             'Last':'Whitney',
             'Age':25,
             'Birthday':'March 10',
             'Grades':AList1,
             'AnArray':np.ndarray.tolist(AnArray),
             'Rvector':np.ndarray.tolist(R),           # Make R into a list to output
             }]



with open("test_dicionary2.json",'w') as fp: 
    variable=json.dumps(vgfin_data, indent=4)
    fp.write(variable)


