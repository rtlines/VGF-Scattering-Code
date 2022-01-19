# functon getkv to read in a k-vector file
#    There are library functions to do this in python alread

# The original code is as follows
# C***********************************************************************
# C*---------------------------------------------------------------------*
#       subroutine getkv(kth,kph,NK,KFILE,ERR,ERRlast,mcount,kcount)    
# C***********************************************************************
# C     subroutine to read in the kvectkors from the file KFILE          *
# C***********************************************************************
# C *** Set the value of NMAX via an included file                     *** 
#        implicit none
#        include 'nmax.inc'
# C****  Variables                 
#        real kth(KNMAX),kph(KNMAX),ERR,ERRlast
#        integer NK, N,mcount,kcount 
#        character*80 KFILE       
# C
#        open(UNIT=45,FILE=KFILE,status='UNKNOWN')
#          read(45,200) NK
#          do N=1,NK                  
#             read(45,210) kth(N),kph(N)
#          end do 
#        read(45,*) ERR, ERRlast, mcount,kcount       
#        close(45)
#        return         
#  200   format(I5)
#  210   format(3f15.7)   
#       end   
#
# This is going to be a pain, because we don't have json to help us
#    We will need some formatted input. 
import numpy as np
#
workfile ="Test_1_10_10.kv"
with open(workfile) as f:
     read_data = f.readline()            # read a whole line
     print(read_data)                    # The first line is just the number of kvectors
     NK=int(read_data)                   # turn the whole line into a number
     print(NK)
     kth = np.zeros(NK)                  # make arrays for the kvector components
     kph = np.zeros(NK)
     khatN = np.zeros((NK,3))              # actual components of the vectors
     for i in range(0,NK,1):             # loop to input all the theta and phi values
          read_data = f.readline()       # read in a whole line
          print(read_data)
          split_string = read_data.split(" ")  #parse the line splitting where
                                               #there are spaces
                                               #Each line starts with a sppace
                                               #then a number, then a space
                                               #then a number. So we will get
                                               #three parts to our split
                                               #a space, a string that is a 
                                               #number, then a second string 
                                               #that is a number.
          print(split_string)
          theta=split_string[1]        # ignore the first split piece and 
                                       # copy the first number string into 
                                       # a variable
          phi=split_string[2]          # and do the same for the second number
                                       # string
          print(theta, phi)     
          kth[i] = float(theta)        # now convert the string into a number
                                       # and put it into the array
          kph[i] = float(phi)          # and do the same for the second string
          print(kth[i],kph[i])
          
     read_data = f.readline()       # The last line has dummy strings
                                    # for kth and kph but has at the end 
     print(read_data)
     split_string = read_data.split(" ")
     print(split_string)
     ERR=int(split_string[7])
     ERRlast=int(split_string[14])
     mcount=int(split_string[25])
     kcount=int(split_string[36])
     print(ERR, ERRlast, mcount, kcount)
     
# Now let's try to modify one of the vectors     
divid = 20
iseed = 1  # not used yet
move1kv(1,khatN,kth,kph,NK,iseed,divd)
     
# now we have the kth and kph values in arrays just like the original
# code wanted.
     
#Now there are fuctions that operate on the kvectors

# The original code is as follows
# C***********************************************************************       
# C*---------------------------------------------------------------------*
#       subroutine move1kv(mcount,khatN,kth,kph,NK,iseed,divd)  
# C***********************************************************************
# C *** Set the value of NMAX via an included file                     *** 
#       implicit none
#       include 'nmax.inc'
# C****  Variables          
#       real khatN(KNMAX,3), kth(KNMAX), kph(KNMAX),PI,PI2,rannum
#       real THinc, PHinc, divd
#       real temp1,temp2
#       integer iseed ,NK, mcount
#       parameter (PI=3.141592654)
#       parameter (PI2=6.28318530718)
#       THinc=PI/divd
#       PHinc=2*THinc      
#       temp1=THinc*rannum(iseed)
#       temp2=PHinc*rannum(iseed)
#       kth(mcount)=kth(mcount)+temp1  !Hinc*rannum(iseed)
#       kph(mcount)=kph(mcount)+temp2   !PHinc*rannum(iseed)
#       if (kth(mcount).gt.PI) then
#          kth(mcount)=2*PI-kth(mcount)
#       end if               
#       if (kph(mcount).gt.PI2) then
#          kph(mcount)=kph(mcount)-PI2
#       end if   
#       if(kth(mcount).lt.0.0) then
#          kth(mcount)=-kth(mcount)
#       end if    
#       if(kph(mcount).lt.0.0) then
#          kph(mcount)=-kph(mcount)
#       end if 
#       call kvectors4(mcount,khatN,kth,kph,NK)
#       return
#       end

# it looks like this changes the k-vector positions randomly for the Monty
#   Carlo part of the code.  So let's look at this in python

# this next part is not tested yet @@@@@@ start here 2022-01-18

def move1kv(mcount,khatN,kth,kph,NK,iseed,divd):
    # mcount seems to be a particular k-vector number
    # khatN seems to be the Cartisian components of the k-vector 
    # kth and kph are our angle arrays we input above
    # NK is the number of k-vectors
    # iseed is a randum number seed (do this later)
    # divid seems to be the amount we are dividing pi by to get our increment
    # start by defining an angular increment abount int the theta direction
    
    THinc = np.pi/divd      #Theta increment
    PHinc = 2*THinc         #Phi direction increment is twice as big
    temmp1 = THinc*random.uniform(-1,1)
    temmp2 = PHinc*random.uniform(-1,1)
    kth[mcount] = kth[mcount]+ temp1
    kph[mcount] = kph[mcount] + temp2
    
    # we have a new k-vector angle, but the change ould make us get strange
    # vectors. We want positive angles, but we also want to adjust the angles
    # either in the positive or the negative direction. This sometimes could 
    # make negagive vector angles. Check for this and correct it if it happens
   
    if (kth[mcount] > np.pi):
        kth[mcount] = 2*np.pi-kth[mcount]
    if (kph[mcount] > np.pi*2):
        kph[mcount] = kph[mcount]-np.pi*2
    if (kth[mcount] < 0.0):
        kth[mcount] = -kth[mcount]
    if (kph[mcount] < 0.0):
        kph[mcount] = -kph[mcount]
    # Now we need to find the components. A function called kvectors4 does that
    # in the original code. Let's just do it here
    # we should worry about whether we are starting arrays at 0 or 1!!!!! @@@@
    khatN[mcount,0] = np.sin(kth[mcount])*np.cos(kph[mcount]) 
    khatN[mcount,1] = np.sin(kth[mcount])*np.sin(kph[mcount])
    khatN[mcount,2] = np.cos(kth[mcount])
    # and we are done. We have a set of k-vectors with one vector randomly
    # moved.
    
    
# function kvectors4 calculates the cartisian components of k-vectors    
    
# The original code is as follows
# C***********************************************************************       
# C*---------------------------------------------------------------------*
#       subroutine kvectors4(N,khatN,kth,kph,NK)                      
# C*---------------------------------------------------------------------* 
# C*    Subroutine to set up the  khatN directions.                      *
# C***********************************************************************       
# C *** Set the value of NMAX via an included file                     *** 
#       implicit none
#       include 'nmax.inc'
# C****  Variables          
#       real khatN(KNMAX,3), kth(KNMAX), kph(KNMAX)
#       integer N ,NK
# C
#       khatN(N,1)=sin(kth(N))*cos(kph(N)) 
#       khatN(N,2)=sin(kth(N))*sin(kph(N))
#       khatN(N,3)=cos(kth(N))
# C
#       return
#       end




# It looks like this isn't used ????
# C***********************************************************************       
# C*---------------------------------------------------------------------*
#       subroutine kvectors3(khatN,kth,kph,NK)                      
# C*---------------------------------------------------------------------* 
# C*    Subroutine to set up the  khatN directions.                      *
# C***********************************************************************       
# C *** Set the value of NMAX via an included file                     *** 
#       implicit none
#       include 'nmax.inc'
# C****  Variables          
#       real khatN(KNMAX,3), kth(KNMAX), kph(KNMAX)
#       integer N, NK
#         do N=1,NK          
#            khatN(N,1)=sin(kth(N))*cos(kph(N)) 
#            khatN(N,2)=sin(kth(N))*sin(kph(N))
#            khatN(N,3)=cos(kth(N))
#         end do
#       return
#       end   
#       