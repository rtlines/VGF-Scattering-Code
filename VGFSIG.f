C***********************************************************************
      program VGFSIG                                                
C***********************************************************************
C**** Program to calculate total cross sections of an arbitrary        *
C****    particle.  The scattering and absorption cross sections are   *
C****    calculated  directly and summed to find the extinction cross  *
C****    section.  The extinction cross section is also calculated     *
C****    using the Optical Theorem, although do to computational error *
C****    the Optical Theorem result may be very wrong!  The program    *
C****    uses the VGF  output files for its input files. Like VGF the  *
C****    assumed time dependence is given by the expression            *
C****    exp(ikr-iwt). Each cross section involves an integration over *
C****    all angles.  Thus the integration can be done in the body     *
C****    frame just as well as the global frame.                       *
C****                                                                  *
C**** Geometry:                                                        *
C****   The particle is aligned with its symmetry axis along           *
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
C**** Author: R. Todd Lines                                            *
C****                                                                  *
C***********************************************************************
C**** Variable diffinitions                                         ****
      implicit none
      include 'nmax.inc'   
C *** loop counters     
      integer mu
      integer I, NDATA, NUSE
C *** Reals    
      real PI, K, K2, K3, DVOL(NMAX), R(NMAX, 3)
      real DEG, alpha, beta, gamma, psi, RAD, W
      real ER, EI, EXR, EXI, EYR, EYI, EZR, EZI, PA2
      real SIGe, SIGa, SIGs, FPIK
      real dsum
      real RRR(3,3),ERR,ERR0
C     complex quantities
      complex E(NMAX, 3), EPS, X, CI
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
      print *, comments
      read(10, *) ERR,ERR0
      print *, ERR,ERR0
      read(10, *) NUSE        
      print *,NUSE
      read(10, *) W, alpha, beta, gamma, psi, RAD 
      print *,W, alpha, beta, gamma, psi, RAD
      print *, 'input cell positions'
C *** read in cell positions and adjusted volume  
C *** read in epsilon and compute chi (X)                            ***
      read(10, *) ER, EI
      print *, ER, EI  
      EPS = cmplx(ER, EI)
      X = (EPS - 1.0)/(4.0*PI)
       dsum=0.0
      do I = 1, NUSE
         read(10, *) R(I, 1), R(I, 2), R(I, 3), DVOL(I)
c         print *, I, R(I, 1), R(I, 2), R(I, 3), DVOL(I)
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
      close(10)
C *** Begin Calculation                                              ***
      write(6, *) 'Calculating ...'
      K = 2.0*PI/W                             
      K2=K*K      
      K3=K2*K
      FPIK=4*PI*K     
C**** Calculate the absorption corss section                        ****
      call sigAbs(SIGa, E, X, DVOL, FPIK, NUSE)
C**** Calculate the extinction cross section using the optical therem **
      call sigExOT(SIGe,K3,K,R,DVOL,X,E,alpha,beta,gamma,psi,NUSE)
C**** Calculate the scattering cross section                        ****
      call simpson (SIGs,R,k,X,DVOL,E,NUSE,RRR)
C *** normalization constant, this seems to be the standard where    ***
C ***   RAD is the symmetry semi axis.                                 *
      PA2 = 1/(PI*RAD*RAD)    
C**** Output the efficiency sig/(pi*a**2) to a file                 ****
C *** open  the output file                                          ***
      open(10, file = OUTFILE, status = 'UNKNOWN')
      write(*,*)
      write(*,*) 'size parameter:    ',K*RAD
      write(*,*) 'Qext (opt. ther):  ',SIGe*PA2
      write(*,*) 'Qext (Qscat+Qabs): ',(SIGs+SIGa)*PA2
      write(*,*) 'Qscat (calcuated): ',SIGs*PA2
      write(*,*) 'Qabs (calculated): ',SIGa*PA2
C
      write(10,*) 'size parameter:    ',K*RAD
      write(10,*) 'Qext (opt. ther):  ',SIGe*PA2
      write(10,*) 'Qext (Qscat+Qabs): ',(SIGs+SIGa)*PA2
      write(10,*) 'Qscat (calcuated): ',SIGs*PA2
      write(10,*) 'Qabs (calculated): ',SIGa*PA2
      close(10)       
 100  format (5f6.4)
      stop
      end
 C***********************************************************************
C*-------------------------------------------------------------------
      subroutine sigAbs(SIGa, E, X, DVOL, FPIK, NUSE)
C*-------------------------------------------------------------------
C***********************************************************************
C*    subroutine to calculate the absorption cross section.            *
C***********************************************************************
C**** Variables                                                     ****
      include 'nmax.inc'
      integer mu, i, NUSE
      real SIGa, EM, DVOL(NMAX), FPIK
      complex E(NMAX,3), X  
C *** Calculate the absorption corss section                         ***
       SIGa=0.0
       do mu=1,NUSE
          EM=0.0
          do i=1,3
             EM=EM+E(mu,i)*conjg(E(mu,i))
          end do      
          SIGa=SIGa+aimag(X)*DVOL(mu)*EM
       end do
       SIGa=FPIK*SIGa     
      return
      end
C***********************************************************************
C*-------------------------------------------------------------------
      subroutine sigExOT(SIGe,K3,K,R,DVOL,X,E,alpha,beta,gamma,psi,NUSE)
C*-------------------------------------------------------------------
C***********************************************************************
C*    Subroutine to calculate the extintion cross section using the    *
C*       Optical theorem.                                              *
C***********************************************************************
C**** Variables                                                     ****
      include 'nmax.inc'
      integer mu,i,NUSE     
      real SIGe,K,K3,R(NMAX,3),DVOL(NMAX)
      real PI,Rhat(3),RRR(3,3),V(3),E0hat(3)
      real alpha,beta,gamma
      complex F(3), temp, CI, X,E(NMAX,3),C
      complex FE(3)
      parameter (PI = 3.141592654, CI = (0.0, 1.0))
C *** Calculate the extinction cross section using the optical therem **
C *    First find the scattering amplitude in the forward direction    *
       call ROT(RRR,alpha,beta,gamma)
C *** Thus E0hat must be in the x-y plane in the lab frame           ***
      V(1)=cos(psi)
      V(2)=sin(psi)
      v(3)=0.0
      call MV(RRR,V,E0hat)
C *** calculate r-hat direction for the given the incident direction ***
       V(1) = 0
       V(2) = 0
       V(3) = 1
       call MV(RRR,V,Rhat)
C ***  Initialize F                                                  ***
       do I = 1, 3
          F(I) = cmplx(0.0, 0.0)
       end do              
C ***  loop over cell number mu                                      ***
       do mu = 1, NUSE                        
C ***     calculate all factors that don't depend on i               ***
          temp = -CI*K*R(mu,3)            
          C = CI*K3*X*cexp(temp)*DVOL(mu)
C ***     calculate the factor that depensd on Emu,i                 ***
          do I=1,2
             FE(I)=E(mu,I)
          end do
          FE(3)=cmplx(0.0,0.0)
C ***     now put the scattering amplitude together                  ***
          do I = 1, 3
             F(I) = F(I)+C*FE(I)
          end do                         
       end do
       SIGe=0.0  
       do I=1,3
          SIGe=SIGe+E0hat(I)*aimag(F(I))
       end do
       SIGe=4*PI*SIGe
      return
      end
C***********************************************************************
C*-------------------------------------------------------------------
      real function FM(TH,PH,R,k,X,DVOL,E,NUSE,RRR)                 
C*-------------------------------------------------------------------
C***********************************************************************
C* Function to do the main calculation of the magnitude of the         *
C*     scattering amplitude. Angles are in radians                     *
C***********************************************************************
C**** Variable diffinitions                                         ****
      implicit none
      include 'nmax.inc'   
C *** loop counters     
      integer mu
      integer I, J, NUSE
C *** Reals    
      real PI, K, K2,K3, DVOL(NMAX), R(NMAX, 3)
      real TH, PH, RDOT, Rhat(3),THhat(3),PHhat(3)
      real RRR(3,3),V(3),HOR,VER     
C     complex quantities
      complex E(NMAX, 3), X, CI, C,temp,fh,fv,TDE,PDE
C *** Parameters
      parameter (PI = 3.141592654, CI = (0.0, 1.0))              
C *** Calculations                                                  ***
c      print *, ' FM '  ,TH,PH
      K2=K*K      
      K3=K2*K
C         print *, 'calculating FM'
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
c         print *, Rhat(1), Rhat(2), Rhat(3)
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
C ***    separate into _H_orizontal and _V_ertical components and    ***
C *        at the same time convert to diff. scat. cross section (/K2) *
         HOR = (Fh*conjg(Fh))/K2
         VER = (Fv*conjg(Fv))/K2     
         FM  = HOR+VER 
c         print *, TH, PH, FM  ,HOR,VER
      return
      end     
   C***********************************************************************
C*-------------------------------------------------------------------
      Subroutine simpson(Res,R,k,X,DVOL,E,NUSE,RRR)                 
C*-------------------------------------------------------------------
C***********************************************************************
C*                                                                     *
C* Subroutine to perform a Simpson's rule integration in two           *
C*    dimensions by performing a one dimensional simpson's rule to     *
C*    all rows of data in a grid, and then performing the Simpson's    *
C*    rule again on the array of results.                              *
C*    Sub Simpson is set up for spherical cordinates theta an phi      *
C*    where the theta is used as the inner most integration.  The      *
C*    integragration is carried out in three sections to try to        *
C*    spend the most time evaluating the integral where it is most     *
C*    important.                                                       *
C***********************************************************************
C**** Variables                                                     ****
      implicit none
      include 'nmax.inc'
C *** Declarations for use in FM
      real R(NMAX, 3),K,DVOL
      complex X(NMAX),E(NMAX, 3)
      integer NUSE
C *** Declaration for use in simpson
      real RRR(3,3)
      real PI
      real TH1,TH2,PH1,PH2
      real Res
      real qsimp1
      parameter (PI = 3.141592654)     
      RES=0.0
C***  Set up limits of integration for each region and integrate.    ***
      write(6,*) 'Integrating region...'
      TH1=0.0                ! 0 degrees
      TH2=pi !30.0*PI/180.0      ! 30 degrees
      PH1=0.0
      PH2=2.0*PI
      RES=qsimp1(TH1,TH2,PH1,PH2,R,K,X,DVOL,E,NUSE,RRR)         
      print *, 'Result so far ', RES
      return
      end
C***********************************************************************
C*-------------------------------------------------------------------
      real function qsimp1(TH1,TH2,PH1,PH2,R,K,X,DVOL,E,NUSE,RRR)   
C*-------------------------------------------------------------------
C***********************************************************************
C *** Declarations for use in FM
      include 'nmax.inc'
      real R(NMAX, 3),K,DVOL
      complex X(NMAX),E(NMAX, 3)
      integer NUSE
C *** Other declarations     
      INTEGER JMAX
      REAL TH1,TH2,PH1,PH2,s,EPS
      real RRR(3,3)
      PARAMETER (EPS=1.e-2, JMAX=200)
CU    USES trapzd1
      INTEGER j
      REAL os,ost,st
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX                
        call trapzd1(TH1,TH2,PH1,PH2,st,j,R,K,X,DVOL,E,NUSE,RRR)
        s=(4.*st-ost)/3.
        if (abs(s-os).lt.EPS*abs(os)) then
           qsimp1=s
           return
        end if
        print *, j,s,os,s-os       
        os=s
        ost=st
11    continue
      pause 'too many steps in qsimp1'
      END
C***********************************************************************
C*-------------------------------------------------------------------
      SUBROUTINE trapzd1(TH1,TH2,PH1,PH2,s,n,R,K,X,DVOL,E,NUSE,RRR)
C*-------------------------------------------------------------------
C***********************************************************************
C *** Declarations for use in FM
      include 'nmax.inc'
      real R(NMAX, 3),K,DVOL
      real RRR(3,3)
      complex X(NMAX),E(NMAX, 3)
      integer NUSE
C *** Other declarations     
      INTEGER n
      REAL TH1,TH2,PH1,PH2,s,qsimp2
      INTEGER it,j
      REAL del,sum,tnm,PH
      if (n.eq.1) then     
        print *, 'call to qsimp2'
        s=0.5*(PH2-PH1)*(qsimp2(TH1,TH2,PH1,R,K,X,DVOL,E,NUSE,RRR)
     &    +qsimp2(TH1,TH2,PH2,R,K,X,DVOL,E,NUSE,RRR))
      else
        it=2**(n-2)
        tnm=it
        del=(PH1-PH2)/tnm
        PH=PH1+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+qsimp2(TH1,TH2,PH,R,K,X,DVOL,E,NUSE,RRR)
          PH=PH+del       
11      continue
        s=0.5*(s+(PH2-PH1)*sum/tnm)  
      endif
      return
      END                                           
C***********************************************************************
C*-------------------------------------------------------------------
      real function qsimp2(TH1,TH2,PH,R,K,X,DVOL,E,NUSE,RRR)
C*-------------------------------------------------------------------
C***********************************************************************
C *** Declarations for use in FM
      include 'nmax.inc'
      real R(NMAX, 3),K,DVOL
      complex X(NMAX),E(NMAX, 3)
      integer NUSE
C *** Other declarations     
      INTEGER JMAX
      REAL TH1,TH2,s,EPS,PH
      real RRR(3,3)
      PARAMETER (EPS=1.e-2, JMAX=200)
CU    USES trapzd2
      INTEGER j
      REAL os,ost,st
c      print *, 'TH ',TH1, TH2,' PH ',ph
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX  
c        print *, 'call trapzd2 ', TH1, TH2,  PH
        call trapzd2(TH1,TH2,st,j,PH,R,K,X,DVOL,E,NUSE,RRR)
        s=(4.*st-ost)/3.
        if (abs(s-os).lt.EPS*abs(os)) then
           qsimp2=s
           return
        end if
        os=s
        ost=st
11    continue
      pause 'too many steps in qsimp2' 
      END
C***********************************************************************
C*-------------------------------------------------------------------
      SUBROUTINE trapzd2(TH1,TH2,s,n,PH,R,K,X,DVOL,E,NUSE,RRR)
C*-------------------------------------------------------------------
C***********************************************************************
C*    Declarations for use in FM
      include 'nmax.inc'
      real R(NMAX, 3),K,DVOL
      complex X(NMAX),E(NMAX, 3)
      integer NUSE
C*    Other declarations
      INTEGER n
      REAL TH1,TH2,s,FM,PH
      real RRR(3,3)
      INTEGER it,j
      REAL del,sum,tnm,TH
      if (n.eq.1) then
        s=0.5*sin(TH)*(TH2-TH1)*(FM(TH1,PH,R,K,X,DVOL,E,NUSE,RRR)
     &    +FM(TH2,PH,R,K,X,DVOL,E,NUSE,RRR))
      else
        it=2**(n-2)
        tnm=it
        del=(TH2-TH1)/tnm
        TH=TH1+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+sin(TH)*FM(TH,PH,R,K,X,DVOL,E,NUSE,RRR)
          TH=TH+del
11      continue
        s=0.5*(s+(TH2-TH1)*sum/tnm)
      endif
      return
      end
C***********************************************************************
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