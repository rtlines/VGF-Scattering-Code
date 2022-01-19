
# Subroutine MV is a matrix vector multipy. Numpy already does this

# The original code is as follows
# C***********************************************************************
# C*---------------------------------------------------------------------*
#       subroutine MV(M,V,U)                                 
# C*---------------------------------------------------------------------*      
# C***********************************************************************
# C*    subroutine to multiply a 3x3 matrix by a vector with three       *
# C*      components                                                     *
# C***********************************************************************
# C****  Variables                                                    ****   
#       implicit none
#        real M(3,3), V(3), U(3), temp 
#        integer i,j      
# C
#        do i=1,3               
#           temp=0.0
#           do j=1,3
#              temp=temp+M(i,j)*V(j)
#           end do
#           U(i)=temp
#        end do
#       return
#       end                                                                  
# C*********************************************************************** 

# We can achieve the same thign with just   

import numpy as np

M = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
V = np.array([1, 2, 3])
U = M.dot(V)
print(U)

# and we don't really need to make it a subroutine or function, we can just
# use numpy where we need the matrix - vector multipy.  We may need to be 
# careful, though, because Fortran didn't distinguish between row and column
# vectors