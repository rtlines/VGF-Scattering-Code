import numpy as np

# teest the multidimensional array notation
#
A=np.zeros((10,3))
for i in range(10):
    for j in range(3):
        A[i][j]=i*j
        
print(A)
        