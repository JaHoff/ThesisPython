# -*- coding: utf-8 -*-
"""
Perturbation analysis secondary file.
Created on Sat Mar  7 17:20:16 2020

@author: mrfan
"""
from os import getcwd
from os import chdir
import sys 
import numpy as np
import matplotlib.pyplot as pp

sys.path.append('../')
chdir('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF
chdir('PerturbationAnalysis/')


#%% 

chdir("../")
cwd = getcwd()
append = cwd + "/Data/"

print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")

# x,y,z,vx,vy,vz, R, V = (np.ndarray([]) for i in range(0,8))

datafolder = cwd + '/Data/PerturbationAnalysis/'

moon = np.ndarray([])

# Load all the data into 3-dimensional arrays
madevars = False

filename = datafolder + "propagationHistory_0_10km.dat"
t, x0, y0, z0, vx0, vy0, vz0, m0  = DI.ImportPropagationHistory(filename,1, True)

filename = datafolder + "propagationHistory_0_10km_nomoon.dat"
t, x1, y1, z1, vx1, vy1, vz1, m1  = DI.ImportPropagationHistory(filename,0, False)

xd, yd, zd, vxd, vyd, vzd = x1-x0, y1-y0, z1-z0, vx1-vx0, vy1-vy0, vz1-vz0


#%% Some figures

PF.Simple3DPlot(x0, y0, z0, title='Orbit with propagated moon', 
                savefig = False, figfolder = '', name = 'placeholder.png', 
                xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                xmod = 1000., ymod = 1000., zmod = 1000., xlim = None, 
                ylim = None , zlim = None,
                fig = None, ax = None)
    
PF.Simple3DPlot(x1, y1, z1, title='Orbit without propagated moon', 
                savefig = False, figfolder = '', name = 'placeholder.png', 
                xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                xmod = 1000., ymod = 1000., zmod = 1000., xlim = None, 
                ylim = None , zlim = None,
                fig = None, ax = None)
    
    
PF.Simple3DPlot(xd, yd, zd, title='Change vector in orbit from planet propagation', 
                savefig = False, figfolder = '', name = 'placeholder.png', 
                xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                xmod = 1000., ymod = 1000., zmod = 1000., xlim = None, 
                ylim = None , zlim = None,
                fig = None, ax = None)