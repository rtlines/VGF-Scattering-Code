
# Subroutine ROT builds a Euler rotation 

# The original code is as follows
# C***********************************************************************
# C*---------------------------------------------------------------------*
#       subroutine ROT(RRR, alpha, beta, gamma)                                 
# C*---------------------------------------------------------------------*      
# C***********************************************************************
# C*    subroutine to return the Euler rotation matrix where the three   *
# C*      rotation angles are alpha, a rotation about the z-axis, beta,  *
# C*      a rotation about the new y-axis, and gamma, a rotation about   *
# C*      the new z-axis.    Not thourghly checked!!!!!!!                *
# C***********************************************************************
# C****  Variables                                                    **** 
#       implicit none
#        real RRR(3,3),alpha, beta, gamma
# C     
#        RRR(1,1)= cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) 
#        RRR(1,2)= sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma)      
#        RRR(1,3)=-sin(beta)*cos(gamma)
# C     
#        RRR(2,1)=-cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma) 
#        RRR(2,2)=-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma) 
#        RRR(2,3)= sin(beta)*sin(gamma)
# C      
#        RRR(3,1)= cos(alpha)*sin(beta)
#        RRR(3,2)= sin(alpha)*sin(beta)
#        RRR(3,3)= cos(beta)
# C
# 
# c      print *, RRR(1,1),RRR(1,2),RRR(1,3)
# c      print *, RRR(2,1),RRR(2,2),RRR(2,3)
# c      print *, RRR(3,1),RRR(3,2),RRR(3,3)
#       return
#       end

# I think we can use the numpy matrix creation for this

import numpy as np

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

alpha = 0 #np.pi/3
beta = np.pi/2
gamma = np.pi/10

RRR = ROT(alpha, beta, gamma)
print(RRR)

