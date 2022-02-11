# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 19:08:43 2022

@author: rtlines
"""

particle_list=[]

from vpython import *
for i in range(3):
    for j in range(3):
        for k in range(3):
            ball = sphere(pos=vector(i,j,k),radius=0.5)
            particle_list.append(ball)
#ball1 = sphere(pos = vector (1,1,1), radius = 0.2)
#ball2 = sphere(pos = vector(1,2,2), radius = 0.5)
particle=compound(particle_list)
particle.axis = vector(0,1,1)
