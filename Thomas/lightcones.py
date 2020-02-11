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
    def plot(self,figure, shift=0, label=None):
        plt.figure(figure)
        plt.plot(self.radius+shift,self.times,label=label)

class Freefall:
    def __init__(self,r0,metric,steps,time,normalized=True):
        self.r0 = r0
        self.metric = metric
        (self.times,self.radius) = metric.fall(r0,steps,time,normalized)
        (self.proper_times,self.proper_radius) = metric.fall_propre(r0,steps,time,normalized)
    def plot(self,figure, shift=0, label=None):
        plt.figure(figure)
        plt.plot(self.radius,self.times+shift,label=label)
    def plot_proper(self,figure, shift=0, label=None):
        plt.figure(figure)
        plt.plot(self.proper_radius,self.proper_times+shift,label=label)

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
    def fall(self,r0,steps,time,normalized):
        step = time / steps
        if not normalized:
            step/=self.Rs
            r0/=self.Rs
        self.r0 = r0
        radius = [r0]
        times = [0]
        integrator = ode(self.fall_eq)
        integrator.set_integrator("vode")
        integrator.set_initial_value(r0-r0*np.finfo(float).epsneg,0)
        for k in range(steps):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def fall_eq(self,t,r):
        return -np.sqrt((1/r-1/self.r0)*(1-1/r)/(1-1/self.r0))
    def fall_propre(self,r0,steps,time,normalized):
        step = time / steps
        if not normalized:
            step/=self.Rs
            r0/=self.Rs
        self.r0 = r0
        radius = [r0]
        times = [0]
        integrator = ode(self.fall_eq_propre)
        integrator.set_initial_value(r0-r0*np.finfo(float).epsneg,0)
        for k in range(steps):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def fall_eq_propre(self,t,r):
        return -np.sqrt((1/r-1/self.r0))
        
        
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

class Lemaitre(Metric):
    def __init__(self,M):
        super().__init__(M)
    def cones(self,cT,step,n,normalized):
        if not normalized:
            step/=self.Rs
            r/=self.Rs
        return






