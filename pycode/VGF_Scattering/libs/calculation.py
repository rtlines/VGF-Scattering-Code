import numpy as np
import random

def kvector_components(theta, phi):
    # Function to find the components of a unit vector given in three dimensions
    #   given the sherical componets angles theta and phi
    #
    # We are going to call this using our kth and kph arrays and put the result
    #    into khatN.  So the call should be something like
    #    khatN(N)=kvectors_components(kth[N], kph[N])
    #
    khat=np.zeros(3)
    khat[0] = np.sin(theta)*np.cos(phi) 
    khat[1] = np.sin(theta)*np.sin(phi)
    khat[2] = np.cos(theta)
    return khat

def kvectors3(kth,kph,NK):                      
    khatN=np.zeros((NK,3))
    for N in range(NK):          
        khatN[N][0]=np.sin(kth[N])*np.cos(kph[N]) 
        khatN[N][1]=np.sin(kth[N])*np.sin(kph[N])
        khatN[N][2]=np.cos(kth[N])
    return khatN

def move1kv(mcount,khatN,kth,kph,NK,divd):
    # mcount is the index of a particular k-vector number
    # khatN is the Cartisian components of the k-vector 
    # kth and kph are our angle arrays we input above
    # NK is the number of k-vectors
    # iseed is a randum number seed (do this later)
    # divid is the amount we are dividing pi by to get our increment
    #
    # start by defining an angular increment abount int the theta direction
    #
    THinc = np.pi/divd                 #Theta increment
    PHinc = 2*THinc                    #Phi direction increment is twice as big
    temp1 = THinc*random.uniform(-1,1) # this is our theta delta
    temp2 = PHinc*random.uniform(-1,1) # this is our phi delta
    #
    kth[mcount] = kth[mcount]+ temp1   #Add the delta theta to the partucuar
                                       #  theta value
    kph[mcount] = kph[mcount] + temp2  #Add the delta phi to the partucuar
                                       #   phi value
    # we have a new k-vector angle, but the change could make us get strange
    # vectors. We want positive angles, but we also want to adjust the angles
    # either in the positive or the negative direction. This sometimes could 
    # make negagive vector angles. Check for this and correct it if it happens
   
    if (kth[mcount] > np.pi):
        kth[mcount] = 2*np.pi-kth[mcount]
    if (kph[mcount] > np.pi*2):
        kph[mcount] = kph[mcount]-np.pi*2
    kth[mcount] = np.abs(kth[mcount])
    kph[mcount] = np.abs(kph[mcount])
    # Now we need to find the components. A function called kvectors4 does that
    # in the original code. kvectors4(mcount, khatN, kth, kph, NK)
    #
    khatN[mcount]=kvector_components(kth[mcount], kph[mcount])
    #
    #or we could just do it here
    # we should worry about whether we are starting arrays at 0 or 1!!!!! @@@@
    #khatN[mcount,0] = np.sin(kth[mcount])*np.cos(kph[mcount]) 
    #khatN[mcount,1] = np.sin(kth[mcount])*np.sin(kph[mcount])
    #khatN[mcount,2] = np.cos(kth[mcount])
    #
    # and we are done. We have a set of k-vectors with one vector randomly
    # moved.
    return kth, kph, khatN 

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

def delta(alpha, beta):
    return int(alpha == beta)

def dd(alpha, i, beta, j):
    return int(alpha == beta and i == j)

def GAM(d,k,EPS):
       GAM=0.0+0.0j
       kd=k*d
       b1=(3./(4.*np.pi))**(2./3.)
       GAM=b1*(kd)**2 + 1j*kd**3/(2.*np.pi)          
       return GAM

def GG(R,k,d,EPS,a,i,b,j):
    # Inputs:
    #    R matrix of dipole locations (Nuse,3)
    #    k the wave number
    #    d is a specific cell weighting (passed in as a single value)
    #    EPS is the complex permitivity
    #    a, i, b, and j are indexes to get the right component of the right 
    #    dopole location.
    # set up the complex temporary variables and zero them
    PHZ = 0.0+0.0j
    t1 = 0.0+0.0j
    t2 = 0.0+0.0j
    temp = 0.0+0.0j
    ci = 0.0+1.0j

    #
    K2=k*k 
    # make a space for dipole displacement vectors 
    Rab = np.zeros(3)
    Rhat = np.zeros(3)
    if b != a:  
        # @@@ Check for error in this branch
        # calculat4e the separation distance between two dipoles
        #   Rmn = Rn - Rm and RMAG = |Rmn|
        Rab[0] = R[a][0] - R[b][0]
        Rab[1] = R[a][1] - R[b][1]
        Rab[2] = R[a][2] - R[b][2]
        # now find the magnitude
        RMAG = Rab[0]**2 + Rab[1]**2 +Rab[2]**2 
        RMAG = RMAG**0.5
        # make a unit vector in the Rmn direction
        Rhat[0] = Rab[0]/RMAG
        Rhat[1] = Rab[1]/RMAG
        Rhat[2] = Rab[2]/RMAG
        # our exponent is i(k dot r) form the exponent
        temp = ci*k*RMAG
        PHZ = np.exp(temp)
        #
        t1 = (K2/RMAG)*(delta(i,j) - Rhat[i]*Rhat[j])
        t2 = (ci*k/RMAG**2 - 1.0/RMAG**3) * (delta(i,j)-3.0*Rhat[i]*Rhat[j])
        G = PHZ * (t1 + t2)
    else:    
        G = 4.0*np.pi*GAM(d,k,EPS)/(3.0*d**3)
    return G

def CPSI(R,khatN,k,mm,N,b) :
    KDR=khatN[N].dot(R[b])
    temp = 1j*mm*k*KDR
    return np.exp(temp)

def Calculate_PHI(An, H, Y, E0, NK, NUSE):
    CPHI=0.0+0.0j
    for N in range(NK):
        for j in range(3):
            for M in range(NK):
                for l in range (3):
                    CPHI = CPHI + np.conjugate(An[M][l]*H[M][l][N][j]*An[N][j])
                    # End l loop
                # End M loop
            # End j loop
        # End N loop
    for j in range (3):
        for N in range (NK):
            CPHI = CPHI - (np.conjugate(Y[N][j])*An[N][j]                     \
                           +Y[N][j]*np.conjugate(An[N][j]))
            # End N loop
        # End j loop
    for a in range(NUSE):
        for i in range(3):
            CPHI - CPHI + E0[a][i]*np.conjugate(E0)[a][i]
            # End i loop
        # End a loop        
    return CPHI

def ecalc(NUSE,R, An,khatN,X,NK,K,mm):
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
    F = np.zeros((NUSE,3),dtype = complex)
    for b in range(0,NUSE):
        for j in range(3):
            F[b][j]=0.0+0.0j
            for N in range(0,NK):
                F[b][j]=F[b][j]+An[N][j]*CPSI(R,khatN,k,mm,N,b)          
    #
    #Done, Calculate the internal E-field to output it.
    #
    for m in range(0, NUSE):
        for I in range(3):
            E[m][I]=F[m][I]/(1.+(4.*np.pi)*X/3.)
    #
    return E

def Calculate_T1_Matrix( D, R, k, khatN, eps, mm, W, NUSE, NK ):
    for a in range(NUSE):
        for i in range(3):
            for N in range (NK):
                  for j in range(3):
                      T1[a][i][N][j] = 0.0+0.0j
                      for b in range (NUSE):
                          dtemp = D[b]
                          T1[a][i][N][j] =  T1[a][i][N][j] +                  \
                                        ( dd(a,i,b,j)-D[b]**3 * W             \
                                        * GG(R,k, dtemp,eps,a,i,b,j))         \
                                        * CPSI(R,khatN, k, mm, N, b)
                          # End of the b loop
                          
                          #print ("T1 loop ",a,i,N,j,T1[a][i][N][j])
                      # End of the j loop         
                  # End of the N loop
            # End of the i loop
        # End of the a loop
    return T1

def Calculate_Y_Matrix(T1,E0, NUSE, NK):
    print('Building the Y matricx')
    for M in range(NK):
          for l in range(3):
                Y[M][l] = 0.0+0.0j    
                for a in range(NUSE):
                      for i in range(3):
                            Y[M][l] = Y[M][l] + np.conjugate(T1[a][i][M][l])  \
                                      *E0[a][i]
                            # print("M, l, a, i",M, l,a,i, T1[a][i][M][l],E0[a][i]  )
                            # End i loop
                        # End a loop
                    # print("M, l, Y[M][l]",M, l, T1[a][i][M][l],E0[a][i]  )
                    # End l loop
    return Y  

def Calculate_H_Matrix(T1, NK):
    print('Building the H matrix')
    for M in range (NK):
        for l in range (3):
            for N in range(NK):
                for j in range(3):
                    H[M][l][N][j] = 0.0+0.0j
                    for a in range(NUSE):
                        for i in range(3):
                            H[M][l][N][j] =  H[M][l][N][j]            \
                                + np.conjugate(T1[a][i][M][l])        \
                                    * (T1[a][i][N][j])
                            # End i loop
                        # End a loop
                    # End j loop
                # End N loop
            # End l looop
        # End M loop
    return H 

def Reform_Y_H_Matricies(Y, H, NK):
    print('reforming the matrix for solving ')
    for n in range(NK):
        for i in range(3):
            m = 3*(n-0)+i     # in fortran arrays start with 1, but 
                              # our python arrays start with 0
            bb[m] = Y[n][i]
            #print('m, i, m, bb ',n, i, m, bb[m])
            for nnp in range (NK):
                    for j in range (3):
                        mp = 3*(nnp-0)+j
                        aa[m][mp] = H[n][i][nnp][j]
                        # print('m, mp, aa[m][mp]',m, mp, aa[m][mp])
                        # End j loop
                    # End nnp loop
            # End i loop
        # End n loop
    return aa, bb


