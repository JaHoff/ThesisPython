# -*- coding: utf-8 -*-
"""
Created on Sun May 24 14:41:01 2020

@author: USER
"""

import numpy as np

# x,y = np.mgrid[0:100,0:100]


# xm  = 0.1
# yf = 10
# zf = 10
# xf,yf,zf = 0.1*x


t = np.arange(0,100,0.01)
x = 0.1 * t
freq = 20
y = np.cos(freq*t/10)
z = np.sin(freq*t*2)


from matplotlib import pyplot as pp 
from mpl_toolkits.mplot3d import Axes3D

fig = pp.figure()
ax = fig.add_subplot(111, projection='3d')
pp.plot(x,y,z)
pp.xlabel('x')
pp.ylabel('y')