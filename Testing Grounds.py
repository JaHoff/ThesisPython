# -*- coding: utf-8 -*-
"""
A quick testing grounds to test the implementation of code in development.

Not annotated, and messy, for obvious reasons
Created on Thu Feb 20 16:15:45 2020

@author: USER
"""

# TESTING GROUND TO SUPPORT DEVELOPMENT

import DataImport as DI
import GeneralFunctions as GF
import numpy as np
import cv2
import gc
import os
import LazyLib as LL
import PlottingFunctions as PF

import time
print("Start synchronizing data files")
# LL.SyncDataFiles()
print("Synchronized data files")
cwd = os.getcwd()
datafolder = cwd + '/Data/TestL4/'
File = datafolder + 'propagationHistory_1yrP.dat'

print("Importing data")
t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,1, True)
print("Imported data succesfully")

#%%
# PF.clearPlots()

# PF.Simple3DPlot(x,y,z, 'Cartesian orbital motion')
# PF.Simple3DPlot(vx,vy,vz, 'Velocity Hodograph', xlabel='Vx [km/s]', ylabel = 'Vy [km/s]', zlabel = 'Vz [km/s]')

sat = np.concatenate((x,y,z), axis=1)
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
baryfig, baryax = PF.Simple3DPlot(co_x,co_y,co_z)
PF.Simple3DScatter(L4b[-1,0],L4b[-1,1],L4b[-1,2], fig=baryfig, ax=baryax)
PF.Simple3DPlot(L4b[:,0], L4b[:,1], L4b[:,2], 'Barycentric frame motion', xlim=(-0,0.4E6), ylim=(-0,0.6E6), fig=baryfig, ax = baryax)

N = co_x.shape[0]
fig2, ax2 = PF.Simple3DPlot(co_x - L4b[:,0].reshape(N,1)*np.ones(co_x.shape), co_y - L4b[:,1].reshape(N,1)*np.ones(co_y.shape)
                , co_z - L4b[:,2].reshape(N,1)*np.ones(co_z.shape))
PF.Simple3DScatter(0,0,0, fig= fig2, ax = ax2, title="Relative motion to L4 in barycentric frame")


cwd = os.getcwd()
#%% 
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=False, highlightfaulty = True, baselineMagnitude = BL_m, 
                   baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)

L4 = GF.ConvertL4Coordinates(moon)
l4_x, l4_y, l4_z = DI.ConvertToRelativeFrame(L4[:,0], L4[:,1], L4[:,2], x,y,z)

PF.Simple3DPlot(l4_x, l4_y, l4_z, title="L4-centralized satellite positions")
#%% 
# normals = GF.sampleSphereGrid(aStep = 6, latStop = 180, lonStop = 179)
# PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], title='Grid distributed sample points along unit sphere')
d = cwd + '/Figures/SampleDistribution/'
LL.CheckDir(d)
normals = GF.sampleSphereFib(360)
PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], 'Fibonnaci lattice sample points along unit sphere',True, d,'Fibbonaci.png')

normals = GF.sampleSphereGrid(np.sqrt(360),-180,180,-180,180)
PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], 'Grid lattice sample points along unit sphere',True, d,'Grid.png')

normals = GF.sampleSphereFib(5)
#%%

figfolder = cwd + '/Figures/uvw PSF/'
refimg = cv2.imread(figfolder + 'SampleSky.jpg',cv2.IMREAD_GRAYSCALE)

I_PSF, I_S = GF.PSFanalysis(normals, BL_c,resolution = 100.)
    
PF.Simple2DImage(I_S, title='Sample function',axes = '100kmc')
PF.Simple2DImage(GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF)), np.uint8), title='Point Spread Function')
gc.collect()
#%% 

## WILL LIKELY NOT SEE CONCRETE DIFFERENCES IN THE PSF CORE UNLESS RESOLUTION IS AROUND A METER!
## USING 100 METER NOW YIELDS LITTLE DIFFERENCES WITH 10 METER IN GENERAL PIXEL SHAPE

xvec = np.array([[1],[0],[0]])
I_PSFx, I_Sx = GF.PSFanalysis(xvec, BL_c, resolution = 100.)
I_PSFx = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFx)), np.uint8)
PF.Simple2DImage(I_Sx, title='Sample function x axis',axes = '100kmc')
PF.Simple2DImage(I_PSFx[950:1050,950:1050], title='PSF function x axis',axes = '100kmc')

I_PSFx2, dump = GF.PSFanalysisSlow(xvec, BL_c,  resolution = 50.)
I_PSFx2 = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFx2)), np.uint8)
PF.Simple2DImage(I_PSFx2[1900:2100,1900:2100], title='PSF x axis high resolution',axes = '100kmc')

#%% 

PF.PlotAroundEarth(x,y,z, moon, L4) # interrupt does not quit when window is closed?

#%% time runtime of functions

t0 = time.time()
GF.PSFanalysis(normals, BL_c,resolution = 100.)
te = time.time() - t0
print("Traditional PSF analysis took %s s" %te)
t0 = time.time()
GF.PSFanalysisSlow(normals, BL_c,resolution = 100.)
te = time.time() - t0
print("Alternate PSF analysis took %s s" %te)

#%%
print(LL.folder_dict["Preface"])
