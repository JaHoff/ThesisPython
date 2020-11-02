# -*- coding: utf-8 -*-
"""
Code for experiments with the surfaces of Hill
Created on Mon Aug 17 13:11:15 2020

@author: mrfan
"""

import numpy as np
import sys 
# import cv2
from os import chdir
from os import getcwd

import LazyLib as LL
import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF
from matplotlib import pyplot as pp

figure_folder = LL.folder_dict["Deployment"]


def C(x,y,z,mu,V):
    C =   Mu(x,y,z,mu) - V**2
    return C

def Mu(x,y,z,mu):    
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2 + z**2)
    mel = x**2 + y**2 + 2*(1-mu)/r1 + 2*mu/r2 
    return mel


lims = np.arange(-1.2,1.2,0.01)
x,y,z = np.meshgrid(lims,lims,lims)

mu = LL.Constants["mu*"]
ER = LL.Constants["Lunar_mean_dist"]
mD = 5.972E24 + 7.34767309E22
tc = np.sqrt(ER**3/(LL.Constants["G"] * mD ))

x0 = 0.5
y0 = 0*np.sqrt(3)/2
z0 = 0
V0 = (50* tc/ER)

c0 = C(x0,y0,z0,mu,V0)
#%%
MuMap = Mu(x,y,z,mu)

#%%
l =int( MuMap.shape[0]/2)
Sel = MuMap[:,:,l]

binary = Sel >= c0

im = binary*255

# fig = pp.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(x[:,:,0],y[:,:,0],Sel)
#%%

# fig = pp.figure
# pp.imshow(im.astype(np.uint8),extent=(-1.2,1.2,-1.2, 1.2))
# pp.title("Accesible area with given initial conditions")

#%% 
CMap = C(x,y,z,mu, V0*np.ones(x.shape))
Sel = CMap[:,:,l]
threshold = 3
binary = Sel <= threshold

im = binary*255

fig = pp.figure
pp.imshow(im.astype(np.uint8),extent=(-1.2,1.2,-1.2, 1.2))
pp.title(f"c <= {threshold}")

#%%
Lx = [0.5,0.5,-1,0.8403997,1.1595926]
Ly = [np.sqrt(3)/2, -np.sqrt(3)/2 , 0, 0, 0]

fig = pp.figure()

Sel[Sel > 10] = 10
pp.contourf(x[:,:,l],y[:,:,l],Sel, levels=[2,3], hatches=['//'], cmap='gray')
pp.scatter(-mu,0,s = 120, label="Earth")
pp.scatter(1-mu,0,s = 60, label = "Moon")

pp.scatter(Lx, Ly, label = "Lagrangian points",color = 'r')
pp.legend()
pp.xlabel('x')
pp.ylabel('y')

#%%

x0 = 0.5
y0 = 0*np.sqrt(3)/2
z0 = 0
V0 = (100* tc/ER)

c0 = C(x0,y0,z0,mu,V0)
Sel = MuMap[:,:,l]

binary = Sel >= c0
fig = pp.figure()
pp.contourf(x[:,:,l],y[:,:,l],Sel, levels=[0,c0], hatches=['//'], cmap='gray')
pp.scatter(-mu,0,s = 120, label="Earth")
pp.scatter(1-mu,0,s = 60, label = "Moon")
pp.scatter(0.5,0, marker='d', label = 'satellite', c = 'k')
pp.scatter(Lx, Ly, label = "Lagrangian points",color = 'r')
pp.legend()
pp.xlabel('x')
pp.ylabel('y')
pp.xlim([-1.2,1.2])
pp.ylim([-1.2,1.2])
pp.savefig(figure_folder + "Hill_surface_50")

#%%
x0 = 0.5
y0 = 0*np.sqrt(3)/2
z0 = 0
V0 = (1050* tc/ER)

c0 = C(x0,y0,z0,mu,V0)
Sel = MuMap[:,:,l]

binary = Sel >= c0
fig = pp.figure()
pp.contourf(x[:,:,l],y[:,:,l],Sel, levels=[0,c0], hatches=['//'], cmap='gray')
pp.scatter(-mu,0,s = 120, label="Earth")
pp.scatter(1-mu,0,s = 60, label = "Moon")
pp.scatter(0.5,0, marker='d', label = 'satellite', c = 'k')
pp.scatter(Lx, Ly, label = "Lagrangian points",color = 'r')
pp.legend()
pp.xlabel('x')
pp.ylabel('y')
pp.xlim([-1.2,1.2])
pp.ylim([-1.2,1.2])
pp.savefig(figure_folder + "Hill_surface_1500")
