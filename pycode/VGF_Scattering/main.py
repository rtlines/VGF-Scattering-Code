#!/bin/env python3
import os

import numpy as np

from libs import calculation
from libs import input
from libs import output
from libs import logs

def main() -> None:
    PATH: str = os.path.abspath(os.path.join(os.path.dirname(__file__), "logs"))
    logger = logs.Logger(PATH)
    # Main function goes here
    #
    # Test reading in the k-vectors
    #
    # array of theta and pi values to read in.  These are the direction angles
    #   for the unit k-vectors.
    #
    # Define Variables
    #
    # Array Limits:  First let's set limits on how big the arrays can be
    KMAX = 2000     # maximum number of k-vectors
    NMAX =2000      # maximum number of dipoles
    MMAX = 2        # maximum number of Monte Carlo itterations (starts at 1?)
    #
    # Convinient consDEG=PI/180.0tants
    DEG=np.pi/180.0     # conversion from degress to radians
    #
    ####### read in the dipole locatons and data from the input file
    #
    # The code needsto input some things.  I will hard code them here for now @@@@@
    ERR0 = 0.05
    divd = 100
    #
    # Now read the input file.  The file name should be an input, but hard code
    #    for now @@@@@
    #
    
    comments = "Not Set"          #is the coment string from vgfin
    NUSE = 0.0              # the number of dipoles to use
    wave = 0.0              # the wavelength in relative uints (as long as all lenght units 
                       #    are the same it doesn't matter, but use microns)
    alpha = 0.0            # The tree Euler angles describing the particle orientation
    beta = 0.0
    gamma = 0.0 
    psi = 0.0                # the polarization angle (I think)
    RAD = 0.0                # is the particle semimajor axis (use microns)
    ER = 0.0                 # the complext components of the particle permitivity
    EI = 0.0
    TD = 0.0                 # the dipole size
    R = np.zeros((NMAX,3))    # a NUSE bcom,wave,alpha,beta,gamma,psi,RAD,ER,EI,TD,R,Dy 3 array of cell locations
    D = np.zeros(NMAX)         # a NUSE array of cell weights
    workfile = "Figure3.in" # @@@ needs to be an input
    logger.log(f"Opening {workfile}")
    NUSE,comments,wave,alpha,beta,gamma,psi,RAD,ER,EI,TD,R,D = input.read_input_file(workfile)
    logger.log(f"Closed {workfile}")
    
    ###### Read in the K-Vectors
    NK = 0                                # number of k-vectors
    ERR = 0                               # Monte Carlo Loop Test Value (MCTV)
    ERRlast = ERR                         # Last value of MCTV to see if we are 
                                          #   converging to a better answer
    kth = np.zeros(KMAX)                  # make arrays for the kvector components
    kph = np.zeros(KMAX)
    # array of k-vector cartisian components.
    khatN = np.zeros((KMAX,3))              # arrat of components of the vectors
    kworkfile ="Figure3.kv"   #@@@ needs to be an input

    logger.log(f"Opening {kworkfile}")
    NK,kth,kph, ERR, ERRlast, mcount, kcount = input.read_kv(kworkfile)       # Function call to get the k-vectors
    logger.log(f"Closed {kworkfile}")
    
    # find the components of the k-vectors
    logger.log(f"Calculating k-vectors...")
    khatN = calculation.kvectors3(kth,kph,NK) 
    logger.log(f"k-vectors calculated")
    ###### Now we have the data input, start calculations
    #
    ## Find the index of refraction, epsilon, and the suseptibility, X
    eps = complex(ER,EI)
    #print(eps)
    mm = np.sqrt(eps)
    #print(mm)
    X = complex(0.0,0.0)
    X = (eps - 1)/(4.0*np.pi)
    #print(X)
    #
    ###### Begin the incident field calculations
    #
    k=2.0*np.pi/wave         # wave number
    #
    print('rotate the input field into particle frame, calc field at each dipole')
    # covert angles from degrees to radians
    alpha = alpha * DEG
    beta = beta * DEG
    gamma = gamma * DEG
    psi = psi * DEG
    #
    # set up the rotation matrix with a call to ROT
    #
    RRR = calculation.ROT(alpha, beta, gamma)
    #print(RRR)
    #
    # K-hat is in the z direction in the lab frame, we need to rotate it into the
    #   particle frame
    V=np.zeros(3)                                 # temporary storage
    V[0] = 0.0                                    
    V[1] = 0.0
    V[2] = 1.0
    khat = RRR.dot(V)
    # Thus E0hat must be in the x-y plane in the lab. Rotate it into the particle
    #   frame     
    V[0] = np.cos(psi)
    V[1] = np.sin(psi)
    V[2] = 0.0
    E0hat = RRR.dot(V)
    print("E0hat ",psi, E0hat)
    #
    # Calculate the W factor    Wcalc=X/(1.+(4.*PI/3.)*X)
    W = X/(1.0+(4.*np.pi/3.)*X)
    #
    ###### Calcuate inital and incident field
    #
    # We need a place to put the calcualted electric field at each dipole location
    #
    E0 = np.zeros((NUSE,3), dtype=complex)
    #
    # we need some complex temporary varuables
    temp = 0+0j
    C = 0+0j
    #
    # first calcualte the r dot khat part of the exponential
    #
    for i in range (NUSE):
        RDK = np.dot(khat, R[i])       # khat dot R
        #print("RDK ", RDK)
        temp = 1j*k*RDK                # i times the wave number times RDK
        #print("temp ",temp)
        C = np.exp(temp)               # take the exponential
        #print ("C ",C)
        E0[i] = C*E0hat                # it is a plane wave, so we need the 
                                       # complex amplitude multiplied by the 
                                       # plane wave exponential
        #print('E0',E0[i],'\n')
    #
    # We are going to calcualte the field in a big loop
    #    Set up the big loop.
    #    We will need spaces for the large matracies to go into
    #
    T1 = np.zeros((NUSE,3,NK,3),dtype = complex)
    Y  = np.zeros((NK,3), dtype = complex)
    H  = np.zeros((NK,3,NK,3), dtype = complex)
    An = np.zeros((NK,3), dtype = complex)
    bb = np.zeros((3*NK), dtype = complex)
    aa = np.zeros((3*NK, 3*NK), dtype = complex)
    xx = np.zeros((3*NK), dtype = complex)
    E = np.zeros((NUSE,3),dtype = complex)

    while (ERR > ERR0) and (mcount < MMAX):

        print('Monte Carlo Loop mcount = ', mcount)
        #ERR came from the k-vector file. We read it in before
        #kcount is the number of kvector tries to loop over. ikcount seems to have 
        #   allowed us to start not at the first k-vector.  I don't see why we would
        #   do that so let's try without it. The kcount loop gives us the number
        #   of random tries for best fit for the field.  But if we get a good one
        #   we drop out of the loop.
        for kcount in range (1):        # should be NK @@@@ put it back after testing T1
            ##################### Calculate the T1 matrix -------------------------
            print('Calculating T1, NK =',NK,' kcount =', kcount)
            T1 = calculation.Calculate_T1_Matrix( D, R, k, khatN, eps, mm, W, NUSE, NK )
            # output the T1 matrix for debugging
            output.print_T1_to_file(T1,NUSE, NK)

            ################# Calculate the Y matrix ------------------------------
            # Calculate the H and Y matricies
            print('Building the Y matricx')
            Y = calculation.Calculate_Y_Matrix(T1,E0, NUSE, NK)
               
            output.print_Y_to_file(Y, NK)
            # Y matrix checked on 02/18/2022
            ################# Calculate the Y matrix ------------------------------
            print('Building the H matrix')
            H = calculation.Calculate_H_Matrix(T1, NK)
            
            output.print_H_to_file(H,NK)
            ################# Now get ready for the matrix inversion---------------
            # Now get ready for the matrix inversion
            # But we need H in a good form for the matrix inverter
            # Try to reform it into a NK*3 by NK*3 matrix
            print('reforming the matrix for solving ')
            [aa, bb] = calculation.Reform_Y_H_Matricies(Y, H, NK)
            
            output.print_bb_to_file(bb, NK)
            output.print_aa_to_file(aa, NK)
            # Now we want to solve the matrix equation aa * xx = bb
            # python is supposed to be able to do this with its linear algebra
            # function linalg.solve(aa,bb)
            print('matrix solve')
            xx = np.linalg.solve(aa,bb)
            # now take our solutino and put it back into matrix component format
            for N in range (NK):
                    for i in range (3):
                        M = 3*(N-0)+i
                        An[N][i] = xx[M]
                        # End i loop
                    # End N lNoop
            # Now we need to test our solution to see how good it is
            # 
            PHI = calculation.Calculate_PHI(An, H, Y, E0, NK, NUSE)
            ERR = PHI * np.conjugate(PHI)
            if (ERR>ERR0):      # Not done with the MC loop, keep going
                   if (ERR > ERRlast):  # Reject: our try isn't so good, reset the k-vectors
                      ERRlast = ERR     # update the ERR history
                      [NK,kth,kph, ERR, ERRlast, mcount, kcount] = input.read_kv(kworkfile)       # Function call to get the k-vectors
                      for N in range(NK):
                          khatN[N]=calculation.kvector_components(kth[N], kph[N])
                   else:               # Acept:  The try was better, keep it
                      ERRlast = ERR    # update the ERR history
                      # because this is an accept, calcualte the fields
                      E = calculation.ecalc(NUSE,R, An,khatN,X,NK,k,mm)
                      #If we want more Monty carlo tries, keep going, change the k-vectors
                      if (kcount < NK):
                          [kth, kph, khatN]=calculation.move1kv(mcount,khatN,kth,kph,NK,divd)
                     
            else:               # Accept the MC try
                     E = calculation.ecalc (NUSE,R, An,khatN,X,NK,k,mm)   
                     break # to exit the kcount loop because we feel we are done
            #
        print('end of the kcount loop')        
        # End of kcount loop                        
        #
        # End the Monte Carlo loop
        mcount = mcount +1  
        print ('End Monte Carlo Loop, mcount =',mcount, 'ERR =',ERR)  
        #
    print ("print out the E-vectors")
    output.print_E_to_file(E, NUSE)            
    # And that is the end of the program
    # Of couse I haven't saved off the E-Fields yet.  So there is not output file
    #  yet.

 

if __name__ == "__main__":
    main()
