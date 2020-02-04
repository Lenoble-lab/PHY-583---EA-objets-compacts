# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:07:01 2020

@author: reill
"""

import matplotlib.pyplot as plt
from lightcones import *
import scipy.constants as cst

M = 10**30
r=[400,800,1475,1500,2000,3000,5000]
step = 5
n=2000

metric = Schwarzschild(M)
test_cones = [LightCone(k,metric,step,n) for k in r]

for cones in test_cones:
    plt.plot(cones.radius,cones.times)

plt.show()
