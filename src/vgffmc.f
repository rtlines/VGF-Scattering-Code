C***********************************************************************
      Program VGF
C***********************************************************************
C**** Program to use the F-field formalizm to solve the scattering     *
C****   problem for a particle of arbitrary shape.  The particle is    *
C****   divided into an array of dipoles on a cubic lattice (by        *
C****   program VGFIN).  The scattering is computed through a plane    *
C****   wave expansion of the field inside the particle.  From this    *
C****   the external field and phase function are calculated (in       *
C****   program VGFPHZ)                                                *
C****                                                                  *
C**** Units:                                                           *
C****   All equations are in Gausian units.                            *
C****   Lengths are all relative, that is, if you input a wavelength   *
C****   of 10 um then all other lengths must be in um. You may use     *
C****   m or cm or furlongs if you wish as long as all lengths are     *
C****   in the same units.                                             *
C****                                                                  *
C**** Time dependence:                                                 *
C****   The VGF code uses the time dependence exp(-iwt).               *
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
C****  Last Modification by Lines 21 APR 97                            *
C****                                                                  * 
C**** List of subroutines and their signatures                         *
C****                                                                  *
c****  subroutine kvectors(NTH,khatN,count)                            *
C****  complex Function GAM(d,k,EPS)                                   *
C****  complex function CPSI(R,KhatN,k,mm,N,b)                         *
C****  complex Function Wcalc(X,d,k,EPS)                               *
C****  subroutine Gcalc(R,k,Gmn,d,EPS,m,n)                             *
C****  real function delta(alpha,beta)                                 *
C****  real function dd(alpha,i,beta,k)                                *
C****  subroutine ROT(RRR, alpha, beta, gamma)                         *
C****  subroutine MV(M,V,U)                                            *
C***********************************************************************
C**** Declarations                                                  ****
      implicit none
C *** Set the value of NMAX via an included file                     ***
      include 'nmax.inc'
C**** Varables 
      integer  ipvt(3*KNMAX)
      complex  z(3*KNMAX)  
      real     rcond
C *** Integers
      integer  NUSE, I, J, mcount, kcount
      integer  N,M,l,a,b,NK,np,mp,iseed,ikcount,itemp
C *** Real
      real     Wave, alpha, beta, gamma, psi, RAD
      real     R(NMAX,3), ER, EI, D(NMAX),  RDK, TD
      real     k, PI, DEG, Khat(3)
      real     RRR(3,3),V(3), E0hat(3)
      real     khatN(KNMAX,3),dtemp
      real     dsum,divd 
      real     kth(KNMAX),kph(KNMAX),ERR,ERRlast,ERR0
C *** Complex
      complex  aa(3*KNMAX,3*KNMAX)
      complex  bb(3*KNMAX)  
      complex  mm, EPS, X, W, C, CI, temp,GAM
      complex  E0(NMAX,3),E(NMAX,3)
      complex T1(NMAX,3,KNMAX,3)
      complex  F(NMAX,3)
      complex  H(KNMAX,3,KNMAX,3),Y(KNMAX,3),An(KNMAX,3)
      complex  PHI
C *** Character
      character*80 INFILE, OUTFILE, COMMENTS,KFILE
C**** Function types
      complex  Wcalc, CPSI, GG, CPHI
      real     dd 
C      real     delta
C**** Typed Parameters ****     
      parameter (PI=3.141592654,DEG=PI/180.0,CI=(0.0,1.0))
C**** DEFINITIONS                                                   ****
      iseed=234564
C *** Read the input and output file names and No. of iterations     ***
      write(*,*)       'Enter input file name:'
      read(*,'(1a80)')  INFILE
      write(*,*)       'Enter output file name:'
      read(*,'(1a80)')  OUTFILE
      write (*,*)       OUTFILE  
      write(*,*)       'Enter input kvector file name:'
      read(*,'(1a80)')  KFILE      
      write(*,*)       'Enter exit criteria ERR0, and increment divisor'
      read(*,*)         ERR0,divd 
C *** Read in the input file                                         ***
      write (*,*)      'Reading input file...'
      open (10,file=INFILE)
      read(10,'(1a80)') Comments   
      read(10, *)       NUSE  
      read(10, *)       Wave,alpha,beta,gamma,psi,RAD                             
C**** CALCULATIONS                                                  ****
      read(10,*)        ER, EI, TD    
      EPS=cmplx(ER,EI)
      mm=sqrt(EPS)
      X=(EPS-1.0)/(4.0*PI)
C**** Input the cell positions
      do I=1,NUSE
         read(10,*) R(I,1),R(I,2),R(I,3),D(I)
         dsum=dsum+D(I)**3 
      end do
      write (*,*) 'total volume = ', dsum
      close(10)
C**** Calculate the K-hat direction and K, the wavenumber in vacuum ****
      write (*,*) 'Calculating incident fields...'
      k=2.0*PI/Wave
      alpha=alpha*DEG
      beta =beta*DEG
      gamma=gamma*DEG                    
      psi  =psi*DEG
C *** set up rotation matrix with call to Rot                        ***
      call Rot(RRR,alpha,beta,gamma)                                    
C *** Khat is in the z direction in the lab frame                    ***      
      V(1)=0.0                                    ! temporary storage
      V(2)=0.0
      V(3)=1.0
      call MV(RRR,V,Khat)                                               
C *** Thus E0hat must be in the x-y plane in the lab frame           ***      
      V(1)=cos(psi)
      V(2)=sin(psi)
      v(3)=0.0
      call MV(RRR,V,E0hat)                                              
C *** Assign incident and initial fields...                          ***
      W=Wcalc(X)                       
      do m = 1,NUSE                              
         RDK =0.0
         do I=1,3
            RDK=RDK+Khat(I)*R(m,I)
         end do
         temp=CI*k*RDK                      
         C=cexp(temp)                       
         do I=1,3
            E0(m,I)=C*E0hat(I)               ! incident E-field    
         end do    
      end do
C**** Setup the khatN directions.
      Write (*,*) 'Define the Kn vectors' 
      
      call getkv(kth,kph,NK,KFILE,ERR,ERRlast,mcount,ikcount) 
      ERR=ERRlast
      call kvectors3(khatN,kth,kph,NK)      
C
C     test subroutines here
      a=1
      i=1
      b=1
      j=1
      dtemp=1.24
      print *, "dtemp,k, EPS, a,i, b,j", dtemp,k, EPS, a, i, b, j
      print *,"GG = ", GG(R,k,dtemp, EPS,a,i,b,j)      
C      print *,"GAM = ", GAM(dtemp,k,EPS) 
C**** HERE WE START THE MONTE CARLO LOOP                             ***
      open(50, file = 'ERRlist', STATUS="UNKNOWN")
      write(50,*) INFILE
      write(50,*) 'ERROR LIST'
      do 10 while ((ERR.gt.ERR0) .and. (mcount.lt.100))
      do 15 kcount=ikcount,NK-1 
C *** Internal field calculation                                     ***
      write(*,*) 'Calculating internal fields...'
      print *, 'calculating T1 '
C      open(unit=43,file='t.dat',status='UNKNOWN')
      do a=1,NUSE
         do i=1,3
            do N =1,NK
C               write(43,22) KhatN(N,1),KhatN(N,2),KhatN(N,3)
               do j=1,3
                  T1(a,i,N,j)=(0.0,0.0)
                  do b=1,NUSE
                     dtemp=d(b) 
                     T1(a,i,N,j)=T1(a,i,N,j)+
     &               (dd(a,i,b,j)-d(b)**3*W*GG(R,k,dtemp,EPS,a,i,b,j))
     &               *CPSI(R,KhatN,k,mm,N,b)
                  end do 
                  print *, a,i,N,j,T1(a,i,N,j)
               end do   
            end do
         end do
      end do
C      close(43)
 22   format(3F17.10)
C
      write (*,*) 'calculate H and Y'
      do 1 M=1,NK
      do 2 l=1,3 
          Y(M,l)=(0.0,0.0) 
          do a=1,NUSE
             do i=1,3
               Y(M,l)=Y(M,l)+conjg(T1(a,i,M,l))*E0(a,i)
             end do      
          end do            
          do 3 N=1,NK
          do 4 j=1,3
          H(M,l,N,j)=(0.0,0.0)
             do a=1,NUSE
                do i=1,3   
                  H(M,l,N,j)=H(M,l,N,j)+conjg(T1(a,i,M,l))*T1(a,i,N,j)
                end do      
             end do 
 4        continue !end do
 3        continue !end do
 2    continue !end do  
 1    continue !end do   
C *** Now do the matrix invertion and solve the system of equaitons
      write(*,*) 'Now solve the set of equations'
      write (*,*) 'OUTFILE',OUTFILE
C *** reform H into a N*3 by N*3 complex matrix *****
      do n=1,NK
         do i=1,3      
            m =3*(n -1)+i         
            bb(m)=Y(n,i)
            do np=1,NK
               do j=1,3
                  mp=3*(np-1)+j
                  aa(m,mp)=H(n,i,np,j)
                  if (aa(m,mp).eq.(0.0,0.0)) then
                      print *, m,mp,aa(m,mp)
                      print *, n,np,i,j,H(n,np,i,j)
                  end if    
               end do
            end do
         end do
      end do       
C ***  
       do i=1,3*KNMAX
          ipvt(i)=0
          z(i)=(0.0,0.0)
       end do           
       rcond=0.0
       call cgeco(aa,3*KNMAX,3*NK,ipvt,rcond,z)
c        call cgeco(a,lda,n,ipvt,rcond,z)
       call cgesl(aa,3*KNMAX,3*NK,ipvt,bb,0)
c        call cgesl(a,lda,n,ipvt,b,job)
      do N=1,NK
          do i=1,3 
             M=3*(N-1)+i
             An(N,i)=bb(M)
          end do   
      end do
C *** Here we end the Monte Carlo Loop by testing the error         ****
      PHI=CPHI(An,H,Y,E0,NK,NUSE) 
c     ERRlast=ERR
      ERR = PHI*conjg(PHI)
      print *, ' test ERR ', ERR,' > ',ERR0
      write (50,*) mcount, ERR     
      if (ERR.gt.ERR0) then
      print *, ' test ERR ', ERR,' > ',ERRlast      
          if (ERR.gt.ERRlast) then
             print *,'reject'
C             call Ecalc(NUSE,R,E,F,An,khatN,X,NK,K,mm)                       
C             call printout(COMMENTS,ERR,ERR0,NUSE,wave,alpha,beta,gamma,
C     &              psi, RAD, EPS,R,D,E,NK,An,mcount,OUTFILE,khatN,divd)
             call getkv(kth,kph,NK,KFILE,ERR,ERRlast,itemp,itemp) 
             ERRlast=ERR
c             do i=1,NK
             call kvectors3(khatN,kth,kph,NK)
c             end do
          else
C            accept  
             print *, 'accept'
             call Ecalc(NUSE,R,E,F,An,khatN,X,NK,K,mm)                       
             call printout(COMMENTS,ERR,ERR0,NUSE,wave,alpha,beta,gamma,
     &              psi, RAD, EPS,R,D,E,NK,An,mcount,OUTFILE,khatN,divd)
             call printkv(kth,kph,NK,KFILE,ERR,ERRlast,mcount,kcount) 
             ERRlast=ERR
          end if
C
          if (kcount.lt.NK-1) then 
             call move1kv(kcount+1,khatN,kth,kph,NK,iseed,divd) 
          end if
      else
          call printkv(kth,kph,NK,KFILE,ERR,ERRlast,mcount,kcount)    
          call Ecalc(NUSE,R,E,F,An,khatN,X,NK,K,mm)
          call printout(COMMENTS,ERR,ERR0,NUSE,wave,alpha,beta,gamma,
     &              psi, RAD, EPS,R,D,E,NK,An,mcount,OUTFILE,khatN,divd) 
          goto 10   ! Agghhh a goto, there must be a better way!
      end if  
      print *, ERR,ERR0, ERRlast ,mcount,kcount
  15  continue ! kcount loop 
      mcount=mcount+1     
  10  continue
C**** DONE WITH THE MONTY CARLO LOOP                                **** 
      close(50)
      print *, ERR, ERR0, mcount  
      stop                         
C      
      end                               
C                                                               
C****** FUNCTIONS AND SUBROUTINES **************************************      
C                                                                      
C***********************************************************************
C*---------------------------------------------------------------------*
      subroutine Ecalc(NUSE,R,E,F,An,khatN,X,NK,K,mm)    
C***********************************************************************
C     subroutine to print out a copy of the kvectkors to the file      *
C       KFILE                                                          *
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
       implicit none
       include 'nmax.inc'
C****  Variables   
      integer NUSE, b,j,m,i,N,NK 
      complex  F(NMAX,3),E(NMAX,3),X,CPSI
      complex  An(KNMAX,3),mm
      real khatN(KNMAX,3),PI,R(NMAX,3),K
       parameter (PI=3.141592654)             
C**** Now form the f's                                              ****
      do b=1,NUSE
         do j=1,3    
            F(b,j)=(0.0,0.0)
            do N=1,NK
               F(b,j)=F(b,j)+An(N,j)*CPSI(R,KhatN,k,mm,N,b)         
            end do   
         end do   
      end do
C**** Done, Calculate the internal E-field to output it.            ****
      do m = 1,NUSE
        do I = 1,3                        
          E(m,I)=F(m,I)/(1.+(4.*PI)*X/3.)
        end do
      end do    
      return
      end
C***********************************************************************
C*---------------------------------------------------------------------*
      subroutine printout(COMMENTS,ERR,ERR0,NUSE,wave,alpha,beta,gamma,
     & psi, RAD, EPS,R,D,E,NK,An,mcount,OUTFILE,khatN,divd)    
C***********************************************************************
C     subroutine to print out a copy of the kvectkors to the file      *
C       KFILE                                                          *
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
       implicit none
       include 'nmax.inc'
C****  Variables 
C *** Integers
      integer  NUSE, J, mcount
      integer NK,m,n
C *** Real
      real     Wave, alpha, beta, gamma, psi, RAD
      real     R(NMAX,3),  D(NMAX),khatN(KNMAX,3)
      real     ERR,ERR0,divd
C *** Complex
      complex  EPS
      complex  E(NMAX,3)
      complex  An(KNMAX,3)
C *** Character
      character*80 OUTFILE, COMMENTS         
C**** And output to file                                            ****
c      write(6, *) 'Writing output file...',OUTFILE
c      write (*,*) 'comments',COMMENTS
      open(10, file = OUTFILE, STATUS="UNKNOWN")
      write(10, *) COMMENTS          
      write(10, *) ERR,ERR0
      write(10, *) NUSE
      write(10, *) wave, alpha, beta, gamma, psi ,RAD    
      write(10, *) real(EPS), aimag(EPS)
      do m = 1, NUSE
        write(10, 100) R(m, 1), R(m, 2), R(m, 3), D(m)**3
      end do
      do m = 1, NUSE
        write(10, 110) (real(E(m,J)), aimag(E(m,J)), J=1,3)
      end do
      do N=1,NK                  
         write(10,105) N, (real(An(N,J)),aimag(An(N,J)), J=1,3)
     &    ,real(An(N,1)*conjg(An(N,1))+An(N,2)*conjg(An(N,2))
     &   +An(N,3)*conjg(An(N,3)))
      end do                  
      do N=1,NK
         write (10,115) khatN(N,1),khatN(N,2),khatN(N,3)
      end do
      write(10,*) 'ERR = ',ERR, 'mcount =', mcount,' divd ',divd
c      print *, OUTFILE, ' written '
      close(10)
      return
 100  format(4(1x, 1g18.8))
 105  format(I5,7F8.3)
 110  format(6(1x, 1g18.8))                             
 115  format(3f17.9)
      end
C***********************************************************************
C*---------------------------------------------------------------------*
      complex function CPHI(An,H,Y,E0,NK,NUSE)
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
       implicit none
       include 'nmax.inc'
C****  Variables          
      complex  E0(NMAX,3)
      complex  H(KNMAX,3,KNMAX,3),Y(KNMAX,3),An(KNMAX,3)
      integer N,j,M,l,a,i,NK,NUSE
C
      CPHI=(0.0,0.0)
      do 5 N=1,NK
      do 6 j=1,3 
      do 7 M=1,NK
      do 8 l=1,3
          CPHI = CPHI+conjg(An(M,l))*H(M,l,N,j)*An(N,j) 
  8   continue 
  7   continue 
  6   continue 
  5   continue
      do j=1,3
         do N=1,NK
            CPHI=CPHI-(conjg(Y(N,j))*An(N,j)+Y(N,j)*conjg(An(N,j)))
         end do
      end do 
      do a=1,NUSE
         do i=1,3
            CPHI=CPHI+E0(a,i)*conjg(E0(a,i))
         end do
      end do
      return
      end      
C***********************************************************************
C*---------------------------------------------------------------------*
      subroutine printkv(kth,kph,NK,KFILE,ERR,ERRlast,mcount,kcount) 
C***********************************************************************
C     subroutine to print out a copy of the kvectkors to the file      *
C       KFILE                                                          *
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
       implicit none
       include 'nmax.inc'
C****  Variables          
       real kth(KNMAX),kph(KNMAX),ERR,ERRlast
       integer NK, N ,mcount,kcount
       character*80 KFILE       
C
       open(UNIT=45,FILE=KFILE,status='UNKNOWN')
         write(45,200) NK
         do N=1,NK                  
            write(45,210) kth(N),kph(N)
         end do        
       write(45,*)  ERR, ERRlast,mcount, kcount
       close(45)
       return         
 200   format(I5)
 210   format(3f15.12)   
      end
C***********************************************************************
C*---------------------------------------------------------------------*
      subroutine getkv(kth,kph,NK,KFILE,ERR,ERRlast,mcount,kcount)    
C***********************************************************************
C     subroutine to read in the kvectkors from the file KFILE          *
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
       implicit none
       include 'nmax.inc'
C****  Variables                 
       real kth(KNMAX),kph(KNMAX),ERR,ERRlast
       integer NK, N,mcount,kcount 
       character*80 KFILE       
C
       open(UNIT=45,FILE=KFILE,status='UNKNOWN')
         read(45,200) NK
         do N=1,NK                  
            read(45,210) kth(N),kph(N)
         end do 
       read(45,*) ERR, ERRlast, mcount,kcount       
       close(45)
       return         
 200   format(I5)
 210   format(3f15.7)   
      end   
C***********************************************************************       
C*---------------------------------------------------------------------*
      subroutine move1kv(mcount,khatN,kth,kph,NK,iseed,divd)  
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
      implicit none
      include 'nmax.inc'
C****  Variables          
      real khatN(KNMAX,3), kth(KNMAX), kph(KNMAX),PI,PI2,rannum
      real THinc, PHinc, divd
      real temp1,temp2
      integer iseed ,NK, mcount
      parameter (PI=3.141592654)
      parameter (PI2=6.28318530718)
      THinc=PI/divd
      PHinc=2*THinc      
      temp1=THinc*rannum(iseed)
      temp2=PHinc*rannum(iseed)
      kth(mcount)=kth(mcount)+temp1  !Hinc*rannum(iseed)
      kph(mcount)=kph(mcount)+temp2   !PHinc*rannum(iseed)
      if (kth(mcount).gt.PI) then
         kth(mcount)=2*PI-kth(mcount)
      end if               
      if (kph(mcount).gt.PI2) then
         kph(mcount)=kph(mcount)-PI2
      end if   
      if(kth(mcount).lt.0.0) then
         kth(mcount)=-kth(mcount)
      end if    
      if(kph(mcount).lt.0.0) then
         kph(mcount)=-kph(mcount)
      end if 
      call kvectors4(mcount,khatN,kth,kph,NK)
      return
      end
C***********************************************************************       
C*---------------------------------------------------------------------*
      real function rannum(iseed)
C*---------------------------------------------------------------------*
C***********************************************************************       
C *** Set the value of NMAX via an included file                     *** 
      implicit none
      include 'nmax.inc'
      real ran1
C**** Variables          
      integer iseed
      rannum=2.0*ran1(iseed)-1.0
      return
      end 
C***********************************************************************       
C*---------------------------------------------------------------------*
      subroutine kvectors4(N,khatN,kth,kph,NK)                      
C*---------------------------------------------------------------------* 
C*    Subroutine to set up the  khatN directions.                      *
C***********************************************************************       
C *** Set the value of NMAX via an included file                     *** 
      implicit none
      include 'nmax.inc'
C****  Variables          
      real khatN(KNMAX,3), kth(KNMAX), kph(KNMAX)
      integer N ,NK
C
      khatN(N,1)=sin(kth(N))*cos(kph(N)) 
      khatN(N,2)=sin(kth(N))*sin(kph(N))
      khatN(N,3)=cos(kth(N))
C
      return
      end
C***********************************************************************       
C*---------------------------------------------------------------------*
      subroutine kvectors3(khatN,kth,kph,NK)                      
C*---------------------------------------------------------------------* 
C*    Subroutine to set up the  khatN directions.                      *
C***********************************************************************       
C *** Set the value of NMAX via an included file                     *** 
      implicit none
      include 'nmax.inc'
C****  Variables          
      real khatN(KNMAX,3), kth(KNMAX), kph(KNMAX)
      integer N, NK
        do N=1,NK          
           khatN(N,1)=sin(kth(N))*cos(kph(N)) 
           khatN(N,2)=sin(kth(N))*sin(kph(N))
           khatN(N,3)=cos(kth(N))
        end do
      return
      end
C***********************************************************************       
C*---------------------------------------------------------------------*
      complex Function GAM(d,k,EPS)                                     
C*---------------------------------------------------------------------*      
C***********************************************************************
C*    Function to calculate the self term contribution termed GAMMA    *
C*      in the IBM write-up. The form for GAM  is taken from the work  *
C*      of B. T. Draine and J. Goodman, Astrophysical Journal, 405:    *
C*      685-697, 1993 March 10. Two other self term calculations are   *
C*      listed here for reference.  In my experience, the Draine and   *
C*      Goodman formulation is the better of the three.                *
C* Goedecke and O'Brian: Note the sign change due the time dependance  *
C*      in IBM being exp(-iwt). This differes from Goedecke and        *
C*      O'Brien's choice.                                              *
C*       GAM=(3./(4.*PI))**(2./3.)*(kd)**2 + CI*kd**3/(2.*PI)          *
C* All The terms in the Goedecke and O'Brien series.  Goedecke and     *
C*      O'Brien expand the exponential in the self term integral and   *
C*      throw away most of the trems.  This is the result if you keep  *
C*      all the terms.                                                 *
C*       real a                                                        *
C*       complex temp                                                  *
C*       a=d*(3./(4.*PI))**(1./3.)                                     *
C*       temp=CI*k*a                                                   *
C*       GAM=2.*((1.-CI*k*a)*cexp(temp)-1.)                            * 
C*  Single value checed 1 Aug 96   formula checked 6 Aug 96            *
C***********************************************************************
C**** Variables                                                     ****  
      implicit none
       complex EPS,CI      
       real k,d,kd
       real PI
       parameter (PI=3.141592654,CI=(0.0,1.0)) 
c       real b1,b2,b3,S
        real b1
       kd=k*d
C *** Drain and Goodman                                              ***
c       b1=-1.8915316
c       b2=0.1648469
c       b3=-1.7700004
c       S=1./5.
c       GAM=(3./(4.*PI))*((b1+EPS*(b2+b3*S))*kd**2+(2.*CI*(kd**3)/3.))
c      b1=0.0
       b1=(3./(4.*PI))**(2./3.)
       GAM=b1*(kd)**2 + CI*kd**3/(2.*PI)          
      return
      end
C***********************************************************************
C*---------------------------------------------------------------------*
      complex function CPSI(R,KhatN,k,mm,N,b)                                 
C*---------------------------------------------------------------------*
C***********************************************************************
C*    Function to calculate the trial funciton expansion functions     *
C*      CPSI=cexp(i*k*khatN.R(b))                                      *
C*      formula checked 6 Aug 96                                       *
C***********************************************************************
C *** Set the value of NMAX via an included file                     *** 
      implicit none
      include 'nmax.inc'
C****  Variables                                                    ****
       real KDR,KhatN(KNMAX,3),R(NMAX,3),k,PI
       complex temp,CI,mm
       integer i,b,N 
       parameter (PI=3.141592654,CI=(0.0,1.0))
C                                   
       KDR=0.0
       do i = 1,3
          KDR=KDR+KhatN(N,i)*R(b,i)
       end do
       temp=CI*mm*k*KDR                   ! random test, shoud be +
       CPSI=cexp(temp)       
       return
       end
C***********************************************************************
C*---------------------------------------------------------------------*
      complex Function Wcalc(X)                                 
C*---------------------------------------------------------------------*
C***********************************************************************
C*    Function to calculate the W-factor                               *
C*      W(nu)=X/(1+(4*pi/3)*X)                                         *
C*      Checked 21 April 97                                            *
C***********************************************************************
C****  Variables                                                    ****   
      implicit none
       complex X,CI
       real PI
       parameter (PI=3.141592654,CI=(0.0,1.0))
C
       Wcalc=X/(1.+(4.*PI/3.)*X)
      return
      end                                                       
C***********************************************************************
C*---------------------------------------------------------------------*
      complex function GG(R,k,d,EPS,a,i,b,j)      
C*---------------------------------------------------------------------*      
C***********************************************************************
C*    Function to calculate the dyadic Green's function for a dipole   *
C*      in the IBM write-up (equation ???)                             *
C*    single value checked 1 Aug 96  Formula Checked 6 Aug 96          *
C***********************************************************************
C *** Set the value of NMAX via an included file                     ***  
      implicit none
      include 'nmax.inc'
C****  Variables                                                    ****
       complex PHZ,t1,t2,temp,CI,GAM,EPS
       real RMAG,Rhat(3),R(NMAX,3),Rab(3),k,K2,d,PI
       real delta
       integer a,b,i,j
       Parameter(PI=3.141592654,CI=(0.0,1.0))
C
       K2=k*k 
C       d3=d**3
       if(b.ne.a) then
C         calculate separation distance Rmn=Rn-Rm and RMAG=|Rmn|
          Rab(1)=R(a,1)-R(b,1)
          Rab(2)=R(a,2)-R(b,2)
          Rab(3)=R(a,3)-R(b,3)
          RMAG=Rab(1)**2+Rab(2)**2+Rab(3)**2    
          RMAG=RMAG**0.5
C         Make a unit vector in the Rmn direction                    ***
          Rhat(1)=Rab(1)/RMAG 
          Rhat(2)=Rab(2)/RMAG
          Rhat(3)=Rab(3)/RMAG
C                          
          temp=CI*k*RMAG  
          PHZ=cexp(temp) 
C         
          t1=(K2/RMAG)*(delta(i,j)-Rhat(i)*Rhat(j)) 
          t2=(ci*k/RMAG**2-1.0/RMAG**3)*(delta(i,j)-3.*Rhat(i)*Rhat(j))
          GG=PHZ*(t1+t2)  
        else
          GG=4.*PI*GAM(d,k,EPS)/(3.0*d**3)
       end if        
       return
      end

C***********************************************************************
C*---------------------------------------------------------------------*
      real function delta(alpha,beta)
C*---------------------------------------------------------------------*      
C***********************************************************************
C     This is q delta function                                         *
C***********************************************************************
C****  Variables              
      implicit none
      integer alpha, beta
      delta=0.0
      if (alpha.eq.beta) delta=1.0
      return
      end

C***********************************************************************
C*---------------------------------------------------------------------*
      real function dd(alpha,i,beta, j)                                       
C*---------------------------------------------------------------------*      
C***********************************************************************
C     This is really two delta functions,                              *
C        delta(alpha,beta) * delta(i,j)                                *
C***********************************************************************
C****  Variables                              
      implicit none
      integer alpha,i, beta,j
      dd=0.0
      if ((alpha.eq.beta).and.(i.eq.j)) dd=1.0
      return
      end
                            
C***********************************************************************
C*---------------------------------------------------------------------*
      subroutine ROT(RRR, alpha, beta, gamma)                                 
C*---------------------------------------------------------------------*      
C***********************************************************************
C*    subroutine to return the Euler rotation matrix where the three   *
C*      rotation angles are alpha, a rotation about the z-axis, beta,  *
C*      a rotation about the new y-axis, and gamma, a rotation about   *
C*      the new z-axis.    Not thourghly checked!!!!!!!                *
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

c      print *, RRR(1,1),RRR(1,2),RRR(1,3)
c      print *, RRR(2,1),RRR(2,2),RRR(2,3)
c      print *, RRR(3,1),RRR(3,2),RRR(3,3)
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
