C***********************************************************************
      program VGFPHZ                                                  
C***********************************************************************
C**** Program to calculate the differential scattering cross section   *
C****    for an arbitrary particle using the electric fields generated *
C****    by the VGF program.  First the scattering amplitude is        *
C****    calculated as a function of the scattering angle theta for a  *
C****    particular value of the azmuthal angle phi.  Theta and phi    *
C****    are lab or global frame coordinates.  The VGF code performs   *
C****    it's calculations in the particles body frame coordinates.    *
C****    Numerical calculation is more efficient if theta and phi are  *
C****    expressed in the body frame cordinates,  thus the subroutine  *
C****    Rot from VGF is emploied to do the transformations. Like VGF  *
C****    the assumed time dependence is given by the expression        *
C****    exp(ikr-iwt).  The program takes as it's input the output     *
C****    files from VGF containing the electric fields.                *
C****                                                                  *
C**** Geometry:                                                        *
C****   The particle is aligned with it's symmetry axis along          *
C****   the Z-axis.  The incident plane wave may be rotated            *
C****   around the particle by the angles alpha, beta, and gamma       *
C****   (corresponding to the Euler angles as given in Arfkin,         *
C****   Mathematical Methods for Physicists, 3nd Ed., Academic Press,  *
C****   NY, 1970). All calculations are done in this 'body frame'      *
C****   The polarization is given relative to the x-axis in the lab or *
C****   global frame. In this frame the incident plane wave traveles   *
C****   in the z-hat direction and the polarization is in the x-y      *
C****   plane. The direction of the electric field is measured by an   *
C****   angle psi from the x-axis.                                     *
C****                                                                  *
C**** Author: Original code by Lines                                   * 
C****                                                                  *
C**** extensive modification by Lines 28 JAN 97                        *
C****  Last Modification by Lines 28 JAN 97 to match theory            *
C***********************************************************************
C**** Variable diffinitions                                         ****
      implicit none
      include 'nmax.inc'    
C *** loop counters      
      integer IANG,mu
      integer I, J, NDATA, NUSE 
C *** Reals     
      real*8 PI, K, K2, K3, DVOL(NMAX), R(NMAX, 3)
      real*8 DEG, alpha, beta, gamma, psi, RDOT, Rhat(3), RAD, W
      real*8 TH, PH, RRR(3,3), V(3),THhat(3),PHhat(3)
      real*8 ER, EI, EXR, EXI, EYR, EYI, EZR, EZI, PA2, HOR, VER
      real*8 dsum,ERR,ERR0
C *** complex quantities
      complex E(NMAX, 3), EPS, X, CI, C,temp
      complex FH, FV ,TDE, PDE
C *** File names and comments
      character*80 INFILE, OUTFILE, COMMENTS
C *** Parameters
      parameter (PI = 3.141592654, DEG = PI/180.0, CI = (0.0, 1.0))
      parameter (NDATA = 181)             
C *** Read input filename from standard in...                        ***
      write(6, *) 'Enter input file name:'
      read(5, '(1a80)') INFILE
      write(6, *) 'Enter output file name:'
      read(5, '(1a80)') OUTFILE
C *** Read input file data... First allert user                      ***
      write(6, *) 'Reading input file...'          
C *** read in file name, comments, wavelength, angles, radius        ***
      open(10, file = INFILE)
      read(10, '(1a80)') COMMENTS 
      read(10, *) ERR,ERR0
      read(10, *) NUSE         
      print *,NUSE
      read(10, *) W, alpha, beta, gamma, psi, RAD  
      print *,W, alpha, beta, gamma, psi, RAD
      print *, 'input cell positions'
      print *,NUSE
C *** read in cell positions and adjusted volume   
C *** read in epsilon and compute chi (X)                            ***

      read(10, *) ER, EI
      EPS = cmplx(ER, EI)
      X = (EPS - 1.0)/(4.0*PI)

      dsum=0.0
      do I = 1, NUSE
         read(10, *) R(I, 1), R(I, 2), R(I, 3), DVOL(I)
         dsum=dsum+DVOL(I)
      end do
      write(*,*) 'total vol =', dsum
      print *, 'input eps',NUSE

      print *, 'input fields',NUSE
C *** read in the field at each cell, three complex components       ***
      do mu = 1, NUSE
         read(10, *) EXR, EXI, EYR, EYI, EZR, EZI
         E(mu, 1) = cmplx(EXR, EXI)
         E(mu, 2) = cmplx(EYR, EYI)
         E(mu, 3) = cmplx(EZR, EZI)
      end do
      close(10)
C *** Calculate rotation matrix                                      ***
      call ROT(RRR, alpha,beta,gamma) 
C *** Calculate the normalized differential cross section            ***
      write(6, *) 'Calculating phase function...'
      K = 2.0*PI/W   
      K2=K*K
      K3=K2*K
C *** normalization constant, this seems to be the standard where    ***
C *      RAD is the symmetry semi axis.                                *
      PA2 = 1/(PI*RAD*RAD)
C *** open  the output file                                          ***
      open(10, file = OUTFILE, status = 'UNKNOWN')
C *** Do the main calculation of the scattering amplitude            ***
      PH=0.0 ! for now
C *** Loop over angle (rad)                                         ***
      do IANG = 0, NDATA-1
         TH = PI*IANG/(NDATA - 1.0)
C ***    define THhat and PHhat, and rotate them into prime frame    ***
         V(1)= cos(TH)*cos(PH)
         V(2)= cos(TH)*sin(PH)
         V(3)=-sin(TH)
         call MV(RRR,V,THhat)
         V(1)=-sin(PH)
         V(2)= cos(PH)
         V(3)= 0     
         call MV(RRR,V,PHhat)
C ***    THhat and PHhat are now really primed,but drop prime        ***         
C ***    calculate r-hat prime direction for the given theta         ***
         V(1) = sin(TH)*cos(PH)
         V(2) = sin(TH)*sin(PH)
         V(3) = cos(TH)     
         call MV(RRR,V,Rhat)
C ***    loop over cell number mu                                    ***
         fh=(0.0,0.0) 
         fv=(0.0,0.0)
         do mu = 1, NUSE                         
            RDOT = 0.0
C ***       calculate rhatprime dot rprime for the exponential       ***
            do j = 1, 3           
               RDOT = RDOT + R(mu, j)*Rhat(j)
            end do                               
C ***       calculate all factors that don't depend on i and j       ***
            temp = -CI*K*RDOT             
            C = CI*K3*X*DVOL(mu)*cexp(temp)
C ***       now put the scattering amplitude together                ***
            TDE=(0.0,0.0)
            PDE=(0.0,0.0)
C ***       calculate THhat dot E and PHhat dot E
            do I = 1, 3               
               TDE=TDE+THhat(I)*E(mu,I)                  
               PDE=PDE+PHhat(i)*E(mu,I)
            end do  
C ***       now put the scattering amplitude together                ***            
            fh=fh+C*TDE
            fv=fv+C*PDE
         end do

c         do I=1,3
c            
c            fv=fv+PHhat(I)*F(I)
c         end do
C ***    separate into _H_orizontal and _V_ertical components and    ***
C *        at the same time convert to diff. scat. cross section (/K2) *
         HOR = (Fh*conjg(Fh))/K2
         VER = (Fv*conjg(Fv))/K2
C ***    normalize and and output to a file                          ***
         write(10, 100) TH/DEG, PH/DEG, VER*PA2, HOR*PA2
100      format(1x, 4(1e12.4, 2x)) 

      end do     
      close(10)
      stop
      end
C***********************************************************************
C*---------------------------------------------------------------------*
      subroutine ROT(RRR, alpha, beta, gamma)                                 
C*---------------------------------------------------------------------*      
C***********************************************************************
C*    subroutine to return the Euler rotation matrix where the three   *
C*      rotation angles are alpha, a rotation about the z-axis, beta,  *
C*      a rotation about the new y-axis, and gamma, a rotation about   *
C*      the new z-axis.                                                *
C***********************************************************************
C**** Veriables                                                     ****
      real*8 RRR(3,3),alpha, beta, gamma
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
C*---------------------------------------------------------------------*
      subroutine MV(M,V,U)                                 
C*---------------------------------------------------------------------*      
C***********************************************************************
C*    subroutine to multiply a 3x3 matrix by a vector with three       *
C*      components                                                     *
C***********************************************************************
C**** Veriables                                                     ****
      real*8 M(3,3), V(3), U(3), temp 
      integer i,j      
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
