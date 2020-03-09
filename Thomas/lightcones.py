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
    def __init__(self,r,metric,steps,time,normalized=True):
        step = time/steps
        self.r=r
        self.metric = metric
        (self.times,self.radius) = metric.cones(r,step,steps,normalized)
    def plot(self,figure, shift=0, label=None):
        plt.figure(figure)
        plt.plot(self.radius+shift,self.times,label=label)

class Freefall:
    def __init__(self,r0,metric,steps,time,proper=False,normalized=True):
        self.r0 = r0
        self.metric = metric
        (self.times,self.radius) = metric.fall(r0,steps,time,normalized)
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
        return -np.sqrt((1/r-1/self.r0)*((1-1/r)**2)/(1-1/self.r0))
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
            if radius[-1] <=0:
                radius[-1]=0
                break
        integrator = ode(self.cone_n)
        integrator.set_initial_value(r,0)
        for k in range(n):            
            radius.insert(0,integrator.integrate(integrator.t+step))
            times.insert(0,integrator.t)
            if radius[0] <=0:
                radius[0]=0
                break
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def cone_p(self,t,r):
        return ((r-1)/(r+1))
    def cone_n(self,t,r):
        return -1
    def fall(self,r0,steps,time,normalized):
        step = time / steps
        if not normalized:
            step/=self.Rs
            r0/=self.Rs
        self.r0 = r0
        radius = [r0]
        times = [0]
        integrator = ode(self.fall_eq)
        integrator.set_initial_value(r0-r0*np.finfo(float).epsneg,0)
        for k in range(steps):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
            if radius[-1] <= 0.00001:
                radius[-1]=0
                break
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def fall_eq(self,t,r):
        return -np.sqrt((1/r-1/self.r0))*(r-1)/(r*np.sqrt(1-1/self.r0)-np.sqrt(1/r-1/self.r0))

class Schwarzschild_properTime(Metric):
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
        return (1-1/r)*np.sqrt(1-1/self.r0)/(1-1/r)
    def cone_n(self,t,r):
        return -(1-1/r)*np.sqrt(1-1/self.r0)/(1-1/r)
    def fall(self,r0,steps,time,normalized):
        step = time / steps
        if not normalized:
            step/=self.Rs
            r0/=self.Rs
        self.r0 = r0
        radius = [r0]
        times = [0]
        integrator = ode(self.fall_eq)
        integrator.set_initial_value(r0-r0*np.finfo(float).epsneg,0)
        for k in range(steps):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
            if radius[-1] <= 0:
                radius[-1]=0
                break
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def fall_eq(self,t,r):
        return -np.sqrt((1/r-1/self.r0))

class Schwarzschild_properTimeMain(Metric):
    def __init__(self,M,ratios):
        super().__init__(M)
        self.ratios = ratios
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
        return (1-1/r)*np.sqrt(1-1/self.r0)/(1-1/r)
    def cone_n(self,t,r):
        return -(1-1/r)*np.sqrt(1-1/self.r0)/(1-1/r)
    def fall(self,r0,steps,time,normalized):
        step = time / steps
        r0=np.array([r0*k for k in self.ratios])
        if not normalized:
            step/=self.Rs
            r0/=self.Rs
        self.r0 = r0
        radius = [r0]
        times = [0]
        integrator = ode(self.fall_eq)
        integrator.set_initial_value(r0-r0*np.finfo(float).epsneg,0)
        for k in range(steps):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
            for k in range(len(self.ratios)):
                if radius[-1][k] <= 0:
                    radius[-1][k]=0
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def fall_eq(self,t,r):
        return np.array([-(self.ratios[k]**(3/2))*np.sqrt(1/r[k]-1/self.r0[k])*((1-(self.ratios[k]**3)/r[k])/(1-1/r[-1]))*np.sqrt((1-1/self.r0[-1])/(1-(self.ratios[k]**3)/self.r0[k])) for k in range(len(self.ratios))])

class Kruskal(Metric):
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
            if radius[-1] <=0:
                radius[-1]=0
                break
        integrator = ode(self.cone_n)
        integrator.set_initial_value(r,0)
        for k in range(n):            
            radius.insert(0,integrator.integrate(integrator.t+step))
            times.insert(0,integrator.t)
            if radius[0] <=0:
                radius[0]=0
                break
        (times,radius)=(np.array(times),np.array(radius))
        if not normalized:
            times = times*self.Rs
            radius = radius*self.Rs
        return (times,radius)
    def cone_p(self,v,u):
        return 1
    def cone_n(self,v,u):
        return -1
    def fall(self,r0,steps,time,normalized):
        step = time / steps
        if not normalized:
            step/=self.Rs
            r0/=self.Rs
        self.r0 = r0
        radius = [r0]
        t0 = np.log(r0-1)
        times = [t0]
        integrator = ode(self.fall_eq)
        integrator.set_initial_value(r0-r0*np.finfo(float).epsneg,t0)
        for k in range(steps):            
            radius.append(integrator.integrate(integrator.t+step))
            times.append(integrator.t)
            if radius[-1] <= 0.00001:
                radius[-1]=0
                break
        print(radius)
        print(times)
        us = []
        vs = []
        for k in range(len(times)):
            us.append((np.exp((times[k]+radius[k])/2)+(radius[k]-1)*np.exp((-times[k]+radius[k])/2))/2)
            vs.append((np.exp((times[k]+radius[k])/2)-(radius[k]-1)*np.exp((-times[k]+radius[k])/2))/2)
#            r = radius[k]
#            t = times[k]            
#            us.append(((r-1)**0.5)*np.exp(r/2)*np.cosh(t/2))
#            vs.append(((r-1)**0.5)*np.exp(r/2)*np.sinh(t/2))
        (us,vs)=(np.array(us),np.array(vs))
        if not normalized:
            vs = vs*self.Rs
            us = us*self.Rs
        return (vs,us)
    def fall_eq(self,t,r):
        return -np.sqrt((1/r-1/self.r0))*(r-1)/(r*np.sqrt(1-1/self.r0)-np.sqrt(1/r-1/self.r0))
    def fall_eq_s(self,t,r):
        return -np.sqrt((1/r-1/self.r0)*((1-1/r)**2)/(1-1/self.r0))

class Lemaitre(Metric):
    def __init__(self,M):
        super().__init__(M)
    def cones(self,cT,step,n,normalized):
        if not normalized:
            step/=self.Rs
            r/=self.Rs
        return






