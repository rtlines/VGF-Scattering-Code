# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 16:41:23 2021

@author: rtlines
"""

from scipy.spatial.transform import Rotation as R
import numpy as np

theta = np.pi/2

M = [[ np.cos(theta), -np.sin(theta)], 
     [np.sin(theta), np.cos(theta)]]

x=[1,0]


print(M)
print(x)
print(np.matmul(M,x))