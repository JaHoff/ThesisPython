# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 13:56:03 2020

Read pop files to display convergence process
@author: USER
"""

import numpy as np
import sys 

from os import chdir
from os import getcwd

from matplotlib import pyplot as pp
from mpl_toolkits.mplot3d import Axes3D

sys.path.append('../')
chdir('../')
import LazyLib as LL
import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF

subfolder = "/20SatsEudoxos/"
figure_folder = LL.folder_dict_funcs["OPT"] + subfolder

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization' + subfolder

#%% Set up data import loop

PF.ClearPlots()

algorithms = ["de1220"]
algorithms_full = ["Differential Evolution"]

n_sats = 20
n_isl = 4
n_lcl = 1
n_pop = 32
sd = 72
n_gen = 50


algo = algorithms[0]
algo_name = algorithms_full[0]

# vectors to store the initial configuration coordinates rel to swarm core:
conf_x, conf_y, conf_z = (np.NaN*np.zeros((6,10)) for i in range(0,3))

init = False

plot_colours = [ 'r', 'g', 'b']

#preload moon position data
moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")


for isl in range(0,n_isl):
    for gen in range(1,n_gen+1):
        

        file_fitness = f"population_{algo}_sd{sd}_sats{n_sats}_nisl{n_isl}_npop{n_pop}_int{n_lcl}_g{gen}_i{isl}.dat" 
        
        fitness = DI.ImportDataFile(datafolder + file_fitness)
        
        if init == False: 
            pop = np.zeros((len(fitness),fitness.shape[1],n_gen,n_isl))
            init=True
            
            
        pop[:n_pop,:,gen-1,isl] = fitness
                
#%%
cores = pop[:,:3,:,:]
velocities = pop[:,3:6,:,:]
sats_x = pop[:,6::3,:,:]
sats_y = pop[:,7::3,:,:]
sats_z = pop[:,8::3,:,:]

markers = ['d', '*', '^', '.']
colours = ['k','r','g', 'b']

#%% Make plot
fig,ax = PF.MakeFig(dimensions='3D')
gen = 0
for i in range(0,n_isl):
    PF.Simple3DScatter(sats_x[0,:,gen,i],sats_y[0,:,gen,i],sats_z[0,:,gen,i], color=colours[i], marker = markers[i], fig=fig,ax=ax)
    PF.Simple3DScatter(cores[0,0,gen,i],cores[0,1,gen,i],cores[0,2,gen,i], color=colours[i], marker = 's', fig=fig,ax=ax)

PF.Simple3DScatter(0,0,0, marker = 'o', color = 'r', fig=fig, ax = ax)

fig,ax = PF.MakeFig(dimensions='3D')
gen = 49
for i in range(0,n_isl):
    PF.Simple3DScatter(sats_x[0,:,gen,i],sats_y[0,:,gen,i],sats_z[0,:,gen,i], color=colours[i], marker = markers[i], fig=fig,ax=ax)
    PF.Simple3DScatter(cores[0,0,gen,i],cores[0,1,gen,i],cores[0,2,gen,i], color=colours[i], marker = 's', fig=fig,ax=ax)
    
PF.Simple3DScatter(0,0,0, marker = 'o', color = 'r', fig=fig, ax = ax)

#%% 

file_history = f"propagationHistory_{algo}_sd{sd}_sats{20}_nisl48_npop{n_pop}_int3_d365_best.dat"
file_core = f"corePosition_{algo}_sd{sd}_sats{5}_nisl48_npop{n_pop}_int3_d365_best.dat"
# Compute initial satellite positions relative to core
t, x, y, z, vx, vy, vz,ddd  = DI.ImportPropagationHistory(datafolder + file_history,0, False)
sat = np.concatenate((x,y,z), axis=1)

#%% 
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
                   baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
#%% 
N_baseline_too_large = np.sum(BL_m > 100e3)
N_baseline_too_small = np.sum(BL_m < 500)
print( f"  this optima has {N_baseline_too_large} instances with too large baselines, \n \
      and {N_baseline_too_small} instances with too small baselines")