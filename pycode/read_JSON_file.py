# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 15:25:27 2021

@author: rtlines
"""
import json
import numpy as np

NSID=10
#R=np.empty((0,3),float)       # Define an empty array called "R" that is N by 3
#D=np.empty(0,float)          # Define an empty array called "D" that is N by 1

with open("test_dicionary3.json",'r') as fp: 
    variable=fp.read()
    data=json.loads(variable)


#print(data)

# Put the data into variables

NUSE=int(data[0]['NUSE'])
w=float(data[0]['Wavelength'])
alpha=float(data[0]['Alpha'])
beta=float(data[0]['Beta'])
gamma=float(data[0]['Gamma'])
psi=float(data[0]['Psi'])
RAD=float(data[0]['RAD'])
EPS=float(data[0]['EPS'])
ER=float(data[0]['ER'])
EI=float(data[0]['EI'])
TD=float(data[0]['TD'])
NSID=int(data[0]['NSID'])
NLSID=int(data[0]['NLSID'])
MR=float(data[0]['MR'])
MI=float(data[0]['MI'])
R=np.array(data[0]['Rarray'])
D=np.array(data[0]['Darray'])
          

print("NUSE", NUSE)
print('Wavelength', w)
print('alpha',alpha)
print('beta',beta)
print('gamma',gamma)
print('psi',psi)
print('RAD',RAD)
print('EPS', EPS)
print('ER',ER)
print('EI',EI)
print('TD',TD)
print('NSID',NSID)
print('NLSID',NLSID)
print('MR',MR)
print('MI',MI)
print('Darray',D)

print("Rarray",R)


