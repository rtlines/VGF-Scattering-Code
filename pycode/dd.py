# function dd is a double Dirac delta function delta(a, b)*delta(i,j) 

# The original code is as follows
# C***********************************************************************
# C*---------------------------------------------------------------------*
#       real function dd(alpha,i,beta, j)                                       
# C*---------------------------------------------------------------------*      
# C***********************************************************************
# C     This is really two delta functions,                              *
# C        delta(alpha,beta) * delta(i,j)                                *
# C***********************************************************************
# C****  Variables                              
#       implicit none
#       integer alpha,i, beta,j
#       dd=0.0
#       if ((alpha.eq.beta).and.(i.eq.j)) dd=1.0
#       return
#       end

def dd(alpha, i, beta, j):
    d=0.0
    if ((alpha == beta) and (i == j)):
        d=1.0
    return d
# function delta is a Dirac delta function

# The original code is as follows
# C***********************************************************************
# C*---------------------------------------------------------------------*
#       real function delta(alpha,beta)
# C*---------------------------------------------------------------------*      
# C***********************************************************************
# C     This is q delta function                                         *
# C***********************************************************************
# C****  Variables              
#       implicit none
#       integer alpha, beta
#       delta=0.0
#       if (alpha.eq.beta) delta=1.0
#       return
#       end
# 
def delta(alpha, beta):
    d=0.0
    if ( alpha == beta ):
        d=1.0
    return d  
      


# test both functions
alpha = 1
beta = 1
i = 1
j = 0

d=dd(alpha, i, beta, j)
print(d)
dta=delta(alpha, beta)
print(dta)