# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 16:39:34 2020

@author: reill
"""

import numpy as np
import matplotlib.pyplot as plt

M = 1e30
G = 6.674e-11
c = 299792458

def P(r,x,R,p):
    value = p*c*c*(np.sqrt(1-2*x*(r/R)**2)-np.sqrt(1-2*x))/(3*np.sqrt(1-2*x)-np.sqrt(1-2*x*(r/R)**2))
    if value>= 0:
        return value
    return None

compacity = np.linspace(0,1,100)

pressure = []
density = []
radius = []
for x in compacity:
    R = G*M/(x*c*c)
    p = 3*M/(4*np.pi*R**3)
    density.append(p)
    pressure.append(P(0,x,R,p))
    radius.append(R)

plt.subplot(311)
plt.plot(compacity,radius)
plt.subplot(312)
plt.plot(compacity,density)
plt.subplot(313)
plt.plot(compacity,pressure)
plt.plot()