# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as pp
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os

cwd = os.getcwd()
figfolder = cwd + '/Figures/Fresnel Diffraction/'

if not os.path.exists(figfolder):
    os.makedirs(figfolder)
    print('made folder')

zst = 1
zstp = 10
Dz  = 0.1
Nz = np.int32( (zstp-zst)/Dz)
z = np.round( np.arange(zst,zstp,Dz) * np.ones((1,Nz)), 2)

fst = 200
fstp = 30e6
Df = 100
Nf = np.int32(np.ceil((fstp-fst)/Df))
F = (np.arange(fst,fstp, Df) * np.ones((1,Nf) )).transpose()

xst = -1
xstp = 1.4
Dx = 0.01
Nx = np.int32(np.ceil((xstp-xst)/Dx))

x = (np.arange(xstp,xst, -Dx) * np.ones((1,Nx) )).transpose()

C = 2.998e8
R = 1737.4e3 #Lunar radius in meters
wl = C/(F*R) #Normalize WL to lunar radius system
omega = - np.sqrt(2/(wl*z))

#%% 

s,c = scipy.special.fresnel(omega)
FresnelSum = 0.5* ( (0.5 + c)**2 + (0.5 + s)**2)

#%% 
f1, z1 = np.meshgrid(z,F)
pp.figure(1)
pp.contour( f1,z1 , FresnelSum, [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1])

#%% Extract the data I actually want to use

Zlevels = [1,2,3,4,5,6,7,8,9]
thresholdlow = .03

threshold = 1e-4

results = np.ones(len(Zlevels))
resultslow = np.ones(len(Zlevels))
j = 0
for v in Zlevels:
    i = np.where(z==v)[1][0]
    results[j] = F[np.where(FresnelSum[:,i] <= threshold)[0][0]]/1e3
    resultslow[j] = F[np.where(FresnelSum[:,i] <= thresholdlow)[0][0]]
    j+=1

#%%
fig, ax = pp.subplots(figsize=(10,3))
moon = pp.Circle((0,0), 1, color='k')
ax.scatter(Zlevels,np.zeros(len(Zlevels))-.1, color = 'r', label='3% threshold frequency [Hz]')
ax.scatter(Zlevels,np.zeros(len(Zlevels))+.1, color = 'b', label='0.01% threshold frequency [kHz]')
ax.set_xlim([0,10])
ax.set_ylim([-1.4, 1.4])
ax.set_xlabel('Distance in lunar radii')
ax.set_ylabel('Distance in lunar radii')

ax.legend(loc='lower right')
ax.add_artist(moon)

for i, txt in enumerate(results):
    ax.annotate(txt, (Zlevels[i], 0), (Zlevels[i], 0.3), color='b')
    
for i, txt in enumerate(resultslow):
    ax.annotate(txt, (Zlevels[i], 0), (Zlevels[i], -0.4), color='r')

pp.savefig(figfolder + 'Diffpattern.png')

#%% Other figure
F = 87.5e3

R = 1737.4e3 #Lunar radius in meters
wl = C/(F*R) #Normalize WL to lunar radius system

omega = x * np.sqrt(2/(wl*z))
s,c = scipy.special.fresnel(omega)
FresnelSum = 0.5* ( (0.5 + c)**2 + (0.5 + s)**2)

Z,X = np.meshgrid(z,x)
fig, ax = pp.subplots(figsize=(10,3))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=.05)
moon = pp.Circle((0,-1), 1, color='w')
# ax.scatter(Zlevels,np.zeros(len(Zlevels))-.1, color = 'r', label='3% threshold frequency [Hz]')
plt = ax.imshow(FresnelSum, vmin=0, vmax=0.01 , extent=(0,10,-1,1.4)
          , cmap = 'nipy_spectral', interpolation='nearest')
fig.colorbar(plt, cax=cax)
ax.set_xlim([0,10])
ax.set_ylim([-1., 1.4])
ax.set_xlabel('Distance in lunar radii')
ax.set_ylabel('Distance in lunar radii')
ax.add_artist(moon)


pp.savefig(figfolder + 'Diff_lower.png')

fig, ax = pp.subplots(figsize=(10,3))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=.05)
moon = pp.Circle((0,-1), 1, color='w')
# ax.scatter(Zlevels,np.zeros(len(Zlevels))-.1, color = 'r', label='3% threshold frequency [Hz]')
plt = ax.imshow(FresnelSum, vmin=0, vmax=1.3 , extent=(0,10,-1,1.4)
          , cmap = 'nipy_spectral')
fig.colorbar(plt, cax=cax)
ax.set_xlim([0,10])
ax.set_ylim([-1., 1.4])
ax.set_xlabel('Distance in lunar radii')
ax.set_ylabel('Distance in lunar radii')
ax.add_artist(moon)
pp.savefig(figfolder + 'Diff_complete.png')