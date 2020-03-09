# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:08:43 2020

@author: reill
"""

import numpy as np
radius=[5,1,0.01]
times=[0,3,20]


k=2

u = ((radius[k]-1)*np.exp((times[k]+radius[k])/2)+np.exp((times[k]-radius[k])/2))/2
v= ((radius[k]-1)*np.exp((times[k]+radius[k])/2)-np.exp((times[k]-radius[k])/2))/2

print(u,v)