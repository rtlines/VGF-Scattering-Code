C***********************************************************************       
C*---------------------------------------------------------------------*
      Program x2kv                        
C*---------------------------------------------------------------------* 
C*    Program to set up a kvector file for the x2 scattering code      *
C*    Author:  Origonal code  by Lines                                 *
C***********************************************************************       
C *** Set the value of NMAX via an included file                     *** 
      implicit none
      include 'nmax.inc'
      character*80 FNAME
C****  Variables          
      real khatN(KNMAX,3),kth(KNMAX),kph(KNMAX),THI(30),TH, PH,PI 
      real delTh, delPh
      integer NTH,NPH,i,j,count ,N
      parameter (PI=3.141592654)   
      write(*,*) 'Enter the number of angles theta (1-5)'
      read (*,*) NTH
      write(*,*) 'Enter the number of angles phi'
      read (*,*) NPH
      write(*,*) ' Enter the file name '
      read(*,'(1a80)') FNAME
        count=1 
        Write (*,*) 'Define the Kn vectors'
        delTh=PI/(NTH+1)
        delPh=2*PI/NPH
c ***   Strait up first       
        kTh(count)=0.0      
        kPH(count)=0.0
        count=count+1
c ***   now the rest        
c        do i=1,NTH         
c           do j=1,NPH
c              kth(count)=THI(i)
c              kPH(count)=j*delPh           
c              print *,count,kTH(count)*180/PI,kPH(count)*180/PI
c              count=count+1
c           end do            
c        end do   
C        
c ***   now the rest        
        do i=1,NTH          
           Th=i*delTH 
           do j=1,NPH
              PH=j*delPH
              kth(count)= Th
              kph(count)= PH                      
              count=count+1
              print *, TH*180/PI,PH*180/PI
           end do
        end do
c ***   now down        
        kTh(count)=PI      
        kPH(count)=0.0
C ***        
        open(UNIT=20,file=FNAME,status='UNKNOWN') 
        write(20,200) count
        do N=1,count
           write (20,210) kth(N),kph(N)
        end do
        close (20)
      stop
 200  format(I5)
 210  format(3f15.7)     
      end
C***********************************************************************    

