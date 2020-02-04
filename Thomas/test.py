# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:07:01 2020

@author: reill
"""

import matplotlib.pyplot as plt
from lightcones import *
import scipy.constants as cst

#M = 10**30
#Rs = 2*cst.G*M/(cst.c**2)
#r=[400,800,1475,1500,2000,3000,5000]
#step = 5
#n=2000

M = 1
r = [0.5,0.8,0.95,1.05,1.2,1.5,2,3,5]
n = 5000
step = 10**-3

plt.figure(1)
metric = Schwarzschild(M)
test_cones = [LightCone(k,metric,step,n) for k in r]

for cones in test_cones:
    plt.plot(cones.radius,cones.times)


plt.figure(2)
metric = EddingtonFinkelstein(M)
test_cones = [LightCone(k,metric,step,n) for k in r]

for cones in test_cones:
    plt.plot(cones.radius,cones.times)

plt.show()
