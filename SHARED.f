C***********************************************************************
C*-------------------------------------------------------------------
      subroutine ROT(RRR, alpha, beta, gamma)                        
C*-------------------------------------------------------------------
C***********************************************************************
C*    subroutine to return the Euler rotation matrix where the three   *
C*      rotation angles are alpha, a rotation about the z-axis, beta,  *
C*      a rotation about the new y-axis, and gamma, a rotation about   *
C*      the new z-axis.                                                *
C***********************************************************************
C****  Variables                                                    ****
       implicit none
       real RRR(3,3),alpha, beta, gamma
C     
       RRR(1,1)= cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) 
       RRR(1,2)= sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma) 
       RRR(1,3)=-sin(beta)*cos(gamma)
C     
       RRR(2,1)=-cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma) 
       RRR(2,2)=-sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma) 
       RRR(2,3)= sin(beta)*sin(gamma)
C      
       RRR(3,1)= cos(alpha)*sin(beta)
       RRR(3,2)= sin(alpha)*sin(beta)
       RRR(3,3)= cos(beta)
C
      return
      end
C***********************************************************************
C*-------------------------------------------------------------------
      subroutine MV(M,V,U)                                 
C*-------------------------------------------------------------------
C***********************************************************************
C*    subroutine to multiply a 3x3 matrix by a vector with three       *
C*      components                                                     *
C***********************************************************************
C****  Variables                                                    ****
       implicit none
       real M(3,3), V(3), U(3), temp 
       integer i,j      
C
       do i=1,3               
          temp=0.0
          do j=1,3
             temp=temp+M(i,j)*V(j)
          end do
          U(i)=temp
       end do
      return
      end                                                              
C***********************************************************************