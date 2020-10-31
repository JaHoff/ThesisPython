# -*- coding: utf-8 -*-
"""
Cost mapping analysis
Created on Tue Jun  2 19:54:30 2020

@author: Jurriez
"""


import numpy as np
import sys 

from os import chdir
from os import getcwd


sys.path.append('../')
chdir('../')
import LazyLib as LL
import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF
import MayaviPlotting as MF

figure_folder = LL.folder_dict["NumericalSim"]


print("Start synchronizing data files")
LL.SyncDataOptimization()
print("Synchronized data files")
cwd = getcwd()
datafolder = cwd + '/Data/costMapping/'

chdir('Cost Mapping')
figfolder = cwd + '/Figures/Cost Mapping/'
LL.CheckDir(figfolder)

#%%
normals = GF.sampleSphereFib(30)

res = 200e3/50

meshgrid = np.mgrid[-100:100:20j, -100:100:20j]
costgrid = np.zeros((20,20))

iterator = 0


do_once = True
basecost =0;
for xc in range(0,20):
    for yc in range(0,20):
        iterator +=1
#        print( f"Now at run {iterator} of 400")
        File = datafolder + f'propagationHistory_5__{xc}_{yc}.dat.dat'
        
#        print("Importing data")
        t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,0, False)
#        print("Imported data succesfully")
    
        if do_once:
            BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x[:,:4],y[:,:4],z[:,:4],
                                                                                            vx[:,:4], vy[:,:4], vz[:,:4], t)
            basecost = np.sum(BL_m > 100e3) + np.sum( BL_m < 500) + np.sum(BLr_m > 1)
            do_once = False
        #%% 
        BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
        
        #I_PSF, I_S,cost = GF.PSFanalysis(normals, BL_c,resolution = res)
        
        #costgrid[xc,yc] = np.sum(cost)
        costgrid[xc,yc] = np.sum(BL_m > 100e3) + np.sum( BL_m < 500) + np.sum(BLr_m > 1) - basecost
        
#%%

fig = MF.MakeFig("Cost function for satellite positions",bgcolor=(1,1,1), fgcolor=(0,0,0))
MF.PlotSurface(fig, meshgrid[0],meshgrid[1],costgrid,cmap = 'copper')
MF.PlotColorSphere(fig, 1, 0,0,np.max(costgrid),R = 3, c = (0,1,0), Nrad = 180, opacity = 0.3)
MF.PlotColorSphere(fig, 1,
                   (x[0,1]-x[0,0])/1e3,
                   (y[0,1]-y[0,0])/1e3,
                   np.max(costgrid),
                   R = 3, c = (0,1,0),
                   Nrad = 180, opacity = 0.3)
MF.AddDefaultAxes(fig,color=(0,0,0), xlabel = 'x [km]', ylabel = 'y [km]', zlabel = 'z [km]')

MF.Title(fig,"Cost for adding a third satellite in x,y domain")

#%%
import matplotlib.pyplot as pp
fig,ax = PF.MakeFig(figsize = (12,4))
cont = pp.contour(meshgrid[0],meshgrid[1],costgrid/1000, 
           [ 6.5, 8,   9.5, 12, 15,20, 25, 28, 30], cmap = 'copper')
pp.clabel(cont  )
pp.scatter(0,0, color= 'k')
pp.scatter(10,0, color='k')
pp.scatter(5,20, color= 'k')
pp.scatter(5,-20, color='k')
#pp.title('Cost for placing a fifth satellite')
pp.xlabel('x [km]')
pp.ylabel('y [km]')
pp.tight_layout()
pp.savefig(figure_folder+'Costcontour_8.png')