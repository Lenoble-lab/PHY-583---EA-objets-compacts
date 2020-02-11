# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:07:01 2020

@author: reill
"""

import matplotlib.pyplot as plt
from lightcones import *
import scipy.constants as cst
import scipy.integrate as sc


#M = 10**30
#Rs = 2*cst.G*M/(cst.c**2)
#r=[400,800,1475,1500,2000,3000,5000]
#step = 5
#n=2000

M = 1
r = [0.5,0.8,0.95,1.05,1.2,1.5,2,3,5]
n = 50000
step = 10**-4
time = 25
steps = 100000
cones = 7

metric = Schwarzschild(M)
test_fall = Freefall(5,metric,steps,time)


plt.figure(1)
plt.plot(test_fall.times,test_fall.radius,label="observateur à l'infini")
#plt.plot(test_fall.proper_times,test_fall.proper_radius,label="temps propre")

#light_cones = [LightCone(test_fall.radius[k*(steps-1)//(cones-1)],metric,step,n) for k in range(cones)]
#
#
#for k in range(cones):
#    cone = light_cones[k]
#    shift = test_fall.times[k*(steps-1)//(cones-1)]
#    plt.plot(cone.times+shift,cone.radius)

plt.ylim(bottom=0)
plt.axhline(y=1, color="black")
plt.title("Chute libre dans une métrique de Schwarzschild")
plt.ylabel("r/Rs")
plt.xlabel("ct/Rs")
plt.legend()


#test_cones = [LightCone(k,metric,step,n) for k in r]

#for cones in test_cones:
#    cones.plot(1,label="r/Rs = " + str(cones.r))


#plt.figure(1)
#plt.xlim(left=0)
#plt.axvline(x=1, color="black")
#plt.title("Chute libre dans une métrique de Schwarzschild")
#plt.xlabel("r/Rs")
#plt.ylabel("ct/Rs")
#plt.legend()

#metric = EddingtonFinkelstein(M)
#test_cones = [LightCone(k,metric,step,n) for k in r]
#
#for cones in test_cones:
#    cones.plot(2,label="r/Rs = " + str(cones.r))
#
#plt.figure(2)
#plt.xlim(left=0)
#plt.axvline(x=1, color="black")
#plt.title("Cones de lumière dans une métrique de Eddington-Finkelstein")
#plt.xlabel("r/Rs")
#plt.ylabel("ct")
#plt.legend()

plt.show()
