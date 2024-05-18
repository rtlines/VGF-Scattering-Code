def print_T1_to_file(T1,NUSE, NK):
    print("Printing T1 data to T1_file.txt")
    with open('T1_file.txt','w') as Tf:
       for a in range(NUSE):
            for i in range(3):
                for N in range (NK):
                    for j in range(3):
                        temp = f"{a}, {i}, {N}, {j}, {T1[a][i][N][j]}\n"
                        Tf.write(temp)            # write a whole line
                        # End of the j loop         
                    # End of the N loop
            # End of the a loopallowed
       # End of file open

def print_Y_to_file(Y, NK):    
    print('Printing Y data to Y_file.txt')
    with open('Y_file.txt','w') as Yf:
        for M in range(NK):
            for l in range(3):
                temp = f"{M}, {l}, {Y[M][l]}\n"
                Yf.write(temp)            # write a whole line
                # End l loop
              # End M loop   
        # Close file

def print_H_to_file(H,NK):      
    print('Printing H data to H_file.txt')
    with open('H_file.txt','w') as Hf:
        for M in range (NK):
            for l in range (3):
                for N in range(NK):
                    for i in range(3):
                        temp = f"{M}, {l}, {N}, {i}, {H[M][l][N][i]}\n"
                        Hf.write(temp)            # write a whole line
                         # End i loop
                     # End N loop
                 # End l looop
             # End M loop
         # Close File    

def print_GG_to_file(R, k, eps):
    print("Printing GG data to GG_file.txt")
    dtemp = 1.24
    with open('GG_file.txt','w') as Gf:
       for a in range(int(NUSE)):
            for i in range(3):
                for b in range (NK):
                    for j in range(3):
                        temp = f"{a}, {i}, {b}, {j}, {GG(R,k,dtemp,eps,a,i,b,j)}\n"
                        Gf.write(temp)            # write a whole line
                        # End of the j loop         
                    # End of the N loop
            # End of the a loop
       # End of file open

def print_bb_to_file(bb, NK):
    print("Printing bb data to bb_file.txt")
    with open('bb_file.txt','w') as bf:
       for i in range(NK*3):
            temp = f"{i}, {bb[i]}\n"
            bf.write(temp)            # write a whole line
            # End of the i loop
       # End of file open

def print_aa_to_file(aa, NK):
    print("Printing aa data to aa_file.txt")
    with open('aa_file.txt','w') as aaf:
       for i in range(NK*3):
           for j in range(NK*3):
               temp = str(i)+", "+str(j)+", "+ str(aa[i][j])+"\n" 
               aaf.write(temp)            # write a whole line
            # End of the i loop
       # End of file open       

def print_E_to_file(E, NUSE):
    print("Printing E vectors at dipoles to E_file.txt")
    with open('E_file.txt','w') as Ef:
        temp = ""
        for i in range(NUSE):
             temp = temp + str(E[i])+"\n"
             Ef.write(temp)
             # End of the i loop
        # End of file open

