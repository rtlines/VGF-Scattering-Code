#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 12:58:38 2022

@author: rtlines
"""

import numpy as np

a = np.zeros((2,2),dtype = complex)
ainv = np.zeros((2,2),dtype = complex)
print(a)
a[0][0] = 2.0+0j
a[0][1] = 4.0+0j
a[1][0] = 3.0+0j
a[1][1] = 5.0+0j

ainv = np.linalg.inv(a)
print(ainv)
b = np.zeros(2,dtype = complex)
b[0] = 16+0j
b[1] = 21+0j
x = np.zeros(2, dtype=complex)
x = np.linalg.solve(a,b)
print(x)


