# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:07:01 2020

@author: reill
"""

import numpy as np
import matplotlib.pyplot as plt
from lightcones import *
import scipy.constants as cst
import scipy.integrate as sc


#M = 10**30
#Rs = 2*cst.G*M/(cst.c**2)
#r=[400,800,1475,1500,2000,3000,5000]
#step = 5
#n=2000

#singularity
x = np.linspace(-10,10,10000)
plt.plot(x,np.sqrt(x**2+1),color='black')
plt.plot(x,-np.sqrt(x**2+1),color='black')




#event horizon
x = np.linspace(-10,10,10000)
plt.plot(x,x,color='black')
plt.plot(x,-x,color='black')





M = 1
r = [0.5,0.8,0.95,1.05,1.2,1.5,2,3,5]

time = 25
steps = 10000
cones = 11

metric = Kruskal(M)
test_fall = Freefall(5,metric,steps,time)



plt.figure(1)
plt.plot(test_fall.radius,test_fall.times)
#plt.plot(test_fall.proper_times,test_fall.proper_radius,label="temps propre")
print("fall done")
steps = len(test_fall.radius)

#
#
r = [1+((k+1)//cones) for k in range(cones)]
i = steps-1
points=[]

for k in r:
    while k > test_fall.radius[i]:
        i-=1
    points.append(i)

print(r)
print(points)

for k in range(cones):
    print(test_fall.radius[k*(steps-1)//(cones)])


light_cones = [LightCone(test_fall.radius[i],metric,steps,time) for i in points]
print(len(light_cones))
    

for k in range(cones):
    cone = light_cones[k]
    shift = test_fall.times[points[k]]
    plt.plot(cone.radius,cone.times+shift)


plt.ylim(bottom=0,top=25)
plt.xlim(left=0,right=10)
plt.axvline(x=1, color="black")
plt.title("Chute libre dans une métrique de Schwarzschild")
plt.xlabel("r/Rs")
plt.ylabel("ct/Rs")
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

#metric = Kruskal(M)
#test_cones = [LightCone(k,metric,steps,time) for k in r]
#
#for cones in test_cones:
#    cones.plot(2,label="r/Rs = " + str(cones.r))

#plt.figure(2)
#plt.xlim(left=0)
#plt.axvline(x=1, color="black")
#plt.title("Cones de lumière dans une métrique de Eddington-Finkelstein")
#plt.xlabel("r/Rs")
#plt.ylabel("ct")
#plt.legend()

plt.show()
