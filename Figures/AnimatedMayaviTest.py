# -*- coding: utf-8 -*-
"""
Testbed for creating the animated orbit plot
Created on Mon Jun  1 14:53:03 2020

@author: Jurriez
"""

import numpy as np
import sys 

from os import chdir
from os import getcwd
sys.path.append('../')
chdir('../')
import LazyLib as LL
# import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF
import MayaviPlotting as MF

#%% 


figure_folder = LL.folder_dict_funcs["MISC"]


print("Start synchronizing data files")
# LL.SyncDataOptimization()
print("Synchronized data files")
cwd = getcwd()
datafolder = cwd + '/Data/Optimization/'
File = datafolder + 'propagationHistory_best.dat'

print("Importing data")
t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,0, False)
moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")
core = DI.ImportSingleColumnFile(datafolder + "corePosition_best.dat")
print("Imported data succesfully")

figfolder = cwd + '/Figures/Optimization/'
LL.CheckDir(figfolder)

chdir('Cost function')

#%%
N = x.shape[0]
# Propagate motion of the core:

delta_p = np.array( (np.mean(np.roll(x,-1,axis=0)-x,axis=1),
                     np.mean(np.roll(y,-1,axis=0)-y,axis=1),
                     np.mean(np.roll(z,-1,axis=0)-z,axis=1) ))
delta_p[:,-1] = np.zeros((3))
delta_p = np.roll(delta_p.T, 1, axis=0)

core_pos = np.zeros((N,3))
core_pos[0,:] = core

for i in range(1,N):
    core_pos[i,:] = core_pos[i-1,:] + delta_p[i,:]

sat = np.concatenate((x,y,z), axis=1)
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)
L4 = GF.ConvertL4Coordinates(moon)

xl4, yl4, zl4 = DI.ConvertToRelativeFrame(L4[:,0], L4[:,1], L4[:,2], x, y, z)

#%% 
MF.PlotAnimatedEarthOrbit( x[:,:1],y[:,:1],z[:,:1], moon, L4, theme="homeworld")

#%%

fig = MF.MakeFig("movement relative to L4 in barycentric frame",theme="homeworld")
MF.Plot3D(fig, core_x-L4b[:,0].reshape(N,1), core_y-L4b[:,1].reshape(N,1), 
          core_z-L4b[:,2].reshape(N,1), lw = 0.2, col = (1,0,0), scale=100e3)

MF.PlotColorSphere(fig, 3)
MF.AddDefaultAxes(fig,xlabel= 'x 100km', ylabel= 'y 100km', zlabel = 'z 100km')
MF.PlotAnimatedParticles(fig,co_x-L4b[:,0].reshape(N,1), co_y-L4b[:,1].reshape(N,1),
                          co_z-L4b[:,2].reshape(N,1), R = 100e3,theme="homeworld", lw=0.2,speed=3)

#%% Movement relative to core
fig = MF.MakeFig("movement relative to swarm core",theme="homeworld")
MF.Plot3D(fig, co_x-core_x, co_y-core_y,
                          co_z-core_z, lw = 0.2, col = (1,0,0), scale=1e3)

MF.PlotColorSphere(fig, 3,c=(1,0,0))
MF.AddDefaultAxes(fig,xlabel= 'x km', ylabel= 'y km', zlabel = 'z km')
MF.PlotAnimatedParticles(fig,co_x-core_x, co_y-core_y,
                          co_z-core_z, R = 1e3,theme="homeworld", lw=0.2,speed=3)