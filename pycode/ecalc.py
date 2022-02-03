# function Ecalc to find the electric field

# The original code is as follows

#       subroutine Ecalc(NUSE,R,E,F,An,khatN,X,NK,K,mm)    
# C***********************************************************************
# C     subroutine to print out a copy of the kvectkors to the file      *
# C       KFILE                                                          *
# C***********************************************************************
# C *** Set the value of NMAX via an included file                     *** 
#        implicit none
#        include 'nmax.inc'
# C****  Variables   
#       integer NUSE, b,j,m,i,N,NK 
#       complex  F(NMAX,3),E(NMAX,3),X,CPSI
#       complex  An(KNMAX,3),mm
#       real khatN(KNMAX,3),PI,R(NMAX,3),K
#        parameter (PI=3.141592654)             
# C**** Now form the f's                                              ****
#       do b=1,NUSE
#          do j=1,3    
#             F(b,j)=(0.0,0.0)
#             do N=1,NK
#                F(b,j)=F(b,j)+An(N,j)*CPSI(R,KhatN,k,mm,N,b)         
#             end do   
#          end do   
#       end do
# C**** Done, Calculate the internal E-field to output it.            ****
#       do m = 1,NUSE
#         do I = 1,3                        
#           E(m,I)=F(m,I)/(1.+(4.*PI)*X/3.)
#         end do
#       end do    
#       return
#       end
#
# New code
#
# Tested:  Not tested yet
#
import numpy as np
#
def ecalc (NUSE,R,E,F,An,khatN,X,NK,K,mm)  :
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
    for b in range(0,NUSE):
        for j in range(3):
            F[b][j]=(0.0,0.0)
            for N in range(0,NK):
                F[b][j]=F[b][j]+An[N][j]*cpsi(R,KhatN,k,mm,N,b)          
    #
    #Done, Calculate the internal E-field to output it.
    #
    for m in range(0, NUSE):
        for I in range(3):
            E[m][I]=F[m][I]/(1.+(4.*np.pi)*X/3.)
    #
    return E
