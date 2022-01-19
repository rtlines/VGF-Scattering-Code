# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 19:09:29 2021

@author: rtlines
"""
a=2.4
b=6.54321
print_string=" test {:5.2f} and test {:5.8f} again"

print(print_string.format(a,b))
print("{:5.2f} {:5.2f}".format(a,b))