# functon rannum gives randum numbers from -1 to 1
#    There are library functions to do this in python alread

# The original code is as follows
# C***********************************************************************       
# C*---------------------------------------------------------------------*
#       real function rannum(iseed)
# C*---------------------------------------------------------------------*
# C***********************************************************************       
# C *** Set the value of NMAX via an included file                     *** 
#       implicit none
#       include 'nmax.inc'
#       real ran1
# C**** Variables          
#       integer iseed
#       rannum=2.0*ran1(iseed)-1.0
#       return
#       end 
#
#
# So to see how to find randum numbers from -1 to 1 we use the random.uniform
#   function in the python random library. Here is an example that calcuates
#   some random numbers and even plots them to show they go from -1 to 1.

import random 
import matplotlib.pyplot as plt
import numpy as np

x=np.zeros(50)
y=np.zeros(50)

for i in range (50):
    x[i] = random.uniform(-1,1)
    
plt.plot(x,y,'ro')
plt.show    