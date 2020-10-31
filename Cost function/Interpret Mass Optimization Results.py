# -*- coding: utf-8 -*-
"""
Created on Sun May 31 14:29:34 2020

@author: mrfan
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

figure_folder = LL.folder_dict_funcs["OPT"]


print("Start synchronizing data files")
LL.SyncDataOptimization()
print("Synchronized data files")
cwd = getcwd()
datafolder = cwd + '/Data/Optimization/'

figfolder = cwd + '/Figures/Optimization/'
LL.CheckDir(figfolder)

algorithms = ["de1220", "sade", "gaco", "pso", "pso_gen"]
algorithms_full = ["Differential Evolution", "Self Adjusting Differential Evolution",
                                            "Particle Swarm Optimization",
                                            "Particle Swarm Optimization Generational"]

n_algos = len(algorithms)
n_gen = 75

for algo in algorithms:
    
    File = datafolder + 'propagationHistory_'+algo+'_best.dat'
    
    print("Importing data")
    t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,0, False)
    moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")
    print("Imported data succesfully")
    
    
    
    init = False
    for i in range(0,n_gen):
        file = datafolder + f'fitness_swarmPropagation_{i}_{i}.dat'
        
        fitness = DI.ImportSingleColumnFile(file)
        
        if init == False: fit = np.zeros((len(fitness),n_gen))
        
        fit[:,i] = fitness
        
        init=True
    
fit_mean = np.mean(fit, axis=0)
fit_max = np.max(fit,axis=0)
fit_min = np.min(fit,axis=0)

gens = np.arange(0,gen)

fig,ax = PF.MakeFig()
PF.Simple2DPlot(gens, fit_mean, xmod = 1., ymod = 1., fig = fig, ax = ax)
PF.Simple2DPlot(gens, fit_min, xmod = 1., ymod = 1., fig = fig, ax = ax)
PF.Simple2DPlot(gens, fit_max, xmod = 1., ymod = 1., fig = fig, ax = ax)
PF.DressFig(fig, ax, title='Optimization progress over generations',
            savefig = False, folder = figfolder, name = 'generational progress.png',
                 xlabel = "Generation", ylabel = "Cost" , ylim = [0,np.max(fit_max)],
                 legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
                 logScale = False, logScaleX = False)


chdir('Cost function')

PF.CloseAll()
#%%
sat = np.concatenate((x,y,z), axis=1)
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
baryfig, baryax = PF.Simple3DPlot(co_x,co_y,co_z)
PF.Simple3DScatter(L4b[-1,0],L4b[-1,1],L4b[-1,2], fig=baryfig, ax=baryax)
PF.Simple3DScatter(LL.Constants["Lunar_mean_dist"]*LL.Constants["mu"],0,
                   0, fig=baryfig, ax=baryax)
PF.Simple3DScatter(0,0,0, fig=baryfig, ax=baryax)
PF.Simple3DPlot(L4b[:,0], L4b[:,1], L4b[:,2], '', xlim=(-0,0.4E6), ylim=(-0,0.6E6), fig=baryfig, ax = baryax)

N = co_x.shape[0]
fig2, ax2 = PF.Simple3DPlot(co_x - L4b[:,0].reshape(N,1)*np.ones(co_x.shape), co_y - L4b[:,1].reshape(N,1)*np.ones(co_y.shape)
                , co_z - L4b[:,2].reshape(N,1)*np.ones(co_z.shape))
PF.Simple3DScatter(0,0,0, fig= fig2, ax = ax2, title="Relative motion to L4 in barycentric frame")

#%% 
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
                   baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
# 

L4 = GF.ConvertL4Coordinates(moon)
l4_x, l4_y, l4_z = DI.ConvertToRelativeFrame(L4[:,0], L4[:,1], L4[:,2], x,y,z)

MF.PlotAroundEarth(x,y,z, moon, L4)
#%% 
N_baseline_too_large = np.sum(BL_m > 100e3)
N_baseline_too_small = np.sum(BL_m < 500)
print( f" This configuration has {N_baseline_too_large} instances with too large baselines, \n \
      and {N_baseline_too_small} instances with too small baselines")
#%% Animated plot around Earth

# MF.PlotAnimatedEarthOrbit( x,y,z, moon, L4, theme="homeworld")

#%% generation improvement plot



#%%
core = DI.ImportSingleColumnFile(datafolder + "corePosition_best.dat")
core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)

co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)

fig = MF.MakeFig("movement relative to swarm core",theme="homeworld")
MF.Plot3D(fig, co_x-core_x, co_y-core_y,
                          co_z-core_z, lw = 0.1, col = (1,0,0), scale=1e3, opacity = 0.5)

MF.PlotColorSphere(fig, 3,c=(1,1,0),opacity = 1)
# MF.AddDefaultAxes(fig,xlabel= 'x km', ylabel= 'y km', zlabel = 'z km')
MF.PlotAnimatedParticles(fig,co_x-core_x, co_y-core_y,
                          co_z-core_z, R = 1e3,theme="default", lw=0.2,speed=3)

#%% First swarm config

fig,ax = PF.MakeFig(None,'3D')
PF.Simple3DScatter(0, 0, 0, fig=fig, ax=ax)
PF.Simple3DScatter(co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0],
                   fig=fig, ax=ax,marker = 'd')
PF.DressFig(fig, ax, 'Intitial swarm configuration relative to swarm core', 
            xlabel = 'km', ylabel = 'km', zlabel = 'km', legendlist = [ "core", "satellites"])

#PF.Animate3DPlotToGIF(fig, ax, "Initial swarm configuration", figfolder)