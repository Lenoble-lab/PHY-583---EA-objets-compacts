# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:57:37 2020

@author: reill
"""

import numpy as np
import matplotlib.pyplot as plt
from lightcones import *
import scipy.constants as cst
import scipy.integrate as sc

M = 10**30
Rs = 2*cst.G*M/(cst.c**2)
r = [0.5,0.8,0.95,1.05,1.2,1.5,2,3,5]

M=1
Rs=1

ratios = [0.25,0.5,0.75,1]

time = 50*Rs
steps = 10000
cones = 11

time_c = 2*Rs
steps_c = 1000

r0=5*Rs
points = 10

falls = Freefall(r0,Schwarzschild_properTimeMain(M, ratios),steps,time)


#for k in range(1,points+1):
#    r = Rs*(k/points)**3
#    i = k-1
#    j=0
#    while falls[i].radius[j] > r and j<len(falls[i].radius)-1:
#        j+=1
#    shift = falls[i].times[j]
#    cone = LightCone(falls[i].radius[j],Schwarzschild(M*(k/points)**3),steps_c,time_c,normalized=False)
#    plt.plot(cone.radius,cone.times+shift,color='grey')

for k in range(len(ratios)):
    
    plt.plot(falls.radius[:,k],falls.times)
#plt.plot(test_fall.proper_times,test_fall.proper_radius,label="temps propre")
print("fall done")
plt.axvline(x=Rs, color="black")

for k in ratios:
    plt.axvline(x=Rs*k**3, color="grey")

plt.show()