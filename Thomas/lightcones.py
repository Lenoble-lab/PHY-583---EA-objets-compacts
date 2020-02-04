# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:46:27 2020

@author: reill
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
from scipy.integrate import ode

class LightCone:
    def __init__(self,r,metric,step,n,normalized=True):
        self.r=r
        self.metric = metric
        (self.times,self.radius) = metric.cones(r,step,n,normalized)

class Metric:
    def __init__(self,M):
        self.M = M
        self.Rs = 2*cst.G*M/(cst.c**2)

class Schwarzschild(Metric):
    def __init__(self,M):
        super().__init__(M)
    def cones(self,r,step,n,normalized):
        if not normalized:
            step/=self.Rs
            r/=self.Rs
        radius = [r]
        times = [0]
        integrator = ode(self.cone_p)
        integrator.set_initial_value(r,0)
        for k in range(n):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
        integrator = ode(self.cone_n)
        integrator.set_initial_value(r,0)
        if r<1 :
            step = -step
        for k in range(n):            
            radius.insert(0,integrator.integrate(integrator.t+step))
            times.insert(0,integrator.t)
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def cone_p(self,t,r):
        return (1-1/r)
    def cone_n(self,t,r):
        return -(1-1/r)
        
class EddingtonFinkelstein(Metric):
    def __init__(self,M):
        super().__init__(M)
    def cones(self,r,step,n,normalized):
        if not normalized:
            step/=self.Rs
            r/=self.Rs
        radius = [r]
        times = [0]
        integrator = ode(self.cone_p)
        integrator.set_initial_value(r,0)
        for k in range(n):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
        integrator = ode(self.cone_n)
        integrator.set_initial_value(r,0)
        for k in range(n):            
            radius.insert(0,integrator.integrate(integrator.t+step))
            times.insert(0,integrator.t)
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def cone_p(self,t,r):
        return ((r-1)/(r+1))
    def cone_n(self,t,r):
        return -1