# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 13:56:03 2020

Read pop files to display convergence process
@author: USER
"""

def vectorProj(v1,v2):
    """Projection of vector v1 onto v2"""
    
    v3 = np.dot(v1,v2)/np.linalg.norm(v2)
    return v3

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

subfolder = "/25satslargespace/"
figure_folder = LL.folder_dict_funcs["OPT"] + subfolder

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization' + subfolder

#%% Set up data import loop

PF.ClearPlots()

threshold = 100

algorithms = ["de1220"]
algorithms_full = ["Differential Evolution"]

n_sats = 25
n_isl = 32
n_lcl = 5
n_pop = 48
sd = 72
gen = 130
n_day = 365


algo = algorithms[0]
algo_name = algorithms_full[0]

init = False

plot_colours = [ 'r', 'g', 'b']

#preload moon position data
moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")
L4 = DI.computeL4Location(moon)


for isl in range(0,n_isl):

        

    file_fitness = f"fitness_{algo}_sd{sd}_sats{n_sats}_nisl{n_isl}_npop{n_pop}_int{n_lcl}_d{n_day}_g{gen}_i{isl}.dat" 
    file_popdata = f"population_{algo}_sd{sd}_sats{n_sats}_nisl{n_isl}_npop{n_pop}_int{n_lcl}_d{n_day}_g{gen}_i{isl}.dat" 
    
    fitness = DI.ImportDataFile(datafolder + file_fitness)
    popdata = DI.ImportDataFile(datafolder + file_popdata)
    
    indexes = np.where(fitness < threshold)[0]
    
    if len(indexes) > 0:
        if init == False:
            pop = np.array(popdata[indexes,:])
            init = True
            
        else:
            for j in indexes:
                pop = np.append(pop, [popdata[j,:]], axis=0)
        
        
        
        
#%%
fig,ax = PF.MakeFig(figsize = (10,10), dimensions = '3D', proj = 'ortho')

PF.Simple3DScatter(pop[:,0], pop[:,1], pop[:,2], fig=fig, ax=ax)
ax.scatter(0,0,0,marker='s', color = 'r')

ax.quiver(pop[:,0]/1e3,pop[:,1]/1e3,pop[:,2]/1e3,pop[:,3],pop[:,4], pop[:,5])

scale = 500
dir_earth = -np.array([L4[0,0],L4[0,1],L4[0,2]])
dir_earth *= scale/np.linalg.norm(dir_earth)
dir_moon = np.array([L4[0,0] - moon[0,0], L4[0,1] - moon[0,1], L4[0,2] - moon[0,2]])
dir_moon *= scale/np.linalg.norm(dir_moon)

ax.quiver(0,0,0,dir_earth[0], dir_earth[1], dir_earth[2],color='r')
ax.quiver(0,0,0,dir_moon[0], dir_moon[1], dir_moon[2], color='y')

ax.set_xlim([-1e3,5e3])
ax.set_ylim([-3e3,3e3])
ax.set_zlim([-3e3,3e3])

#%% Convert results to barycentric frame
c = popdata.shape[0]

R = np.linalg.norm(moon[0,:3])

moon_t0 = [moon[0,:]]
L4_t0 = [L4[0,:]]
for i in range(1,c):
    moon_t0 = np.append(moon_t0, [moon[0,:]], axis=0)
    L4_t0 = np.append(L4_t0, [L4[0,:]], axis=0)

b_pos = DI.ConvertToBarycentricFrame(moon_t0, popdata[:,:3]+L4_t0)
b_velx, b_vely, b_velz = np.empty_like(b_pos[0]),np.empty_like(b_pos[0]),np.empty_like(b_pos[0])

d,d,d,l4b = DI.ConvertToBarycentricFrame(moon_t0, popdata[:,:3]+L4_t0, exportL4 = True)

# manual velocity conversion

b_x = moon[0,:3]/np.linalg.norm(moon[0,:3])
b_z = np.cross(b_x,moon[0,3:6]/np.linalg.norm(moon[0,3:6]))
b_y = np.cross(b_z,b_x)

for j in range(0,c):
    b_velx[j] = vectorProj(popdata[j,3:6],b_x)
    b_vely[j] = vectorProj(popdata[j,3:6],b_y)
    b_velz[j] = vectorProj(popdata[j,3:6],b_z)

b_vel = np.concatenate([b_velx, b_vely, b_velz], axis=1)
b_vel /= np.linalg.norm(b_vel,axis=1).reshape(c,1)

b_vel /= 1000

fig,ax = PF.MakeFig(figsize = (10,10), dimensions = '3D', proj = 'ortho')

PF.Simple3DScatter(b_pos[0][:]/R, b_pos[1]/R, b_pos[2]/R, fig=fig, ax=ax, xmod =1, ymod = 1, zmod = 1)
#ax.scatter(l4b[0,0],l4b[0,1],l4b[0,2],marker='s', color = 'r')

ax.quiver(b_pos[0][:,0]/R,b_pos[1][:,0]/R,b_pos[2][:,0]/R,b_vel[:,0], b_vel[:,1], b_vel[:,2])

#ax.scatter(LL.Constants["mu*"],0,0, marker='*')
#ax.scatter(1-LL.Constants["mu*"],0,0, marker='*')
ax.scatter(0.5-LL.Constants["mu*"],np.sqrt(3)/2,0, marker = 's', color = 'r')

pp.title("barycentric initial core positions and velocities")
ax.set_zlim([-.0051,.0051])

#%%


    
# #%%
# cores = pop[:,:3,:,:]
# velocities = pop[:,3:6,:,:]
# sats_x = pop[:,6::3,:,:]
# sats_y = pop[:,7::3,:,:]
# sats_z = pop[:,8::3,:,:]

# markers = ['d', '*', '^', '.']
# colours = ['k','r','g', 'b']

# #%% Make plot
# fig,ax = PF.MakeFig(dimensions='3D')
# gen = 0
# for i in range(0,n_isl):
#     PF.Simple3DScatter(sats_x[0,:,gen,i],sats_y[0,:,gen,i],sats_z[0,:,gen,i], color=colours[i], marker = markers[i], fig=fig,ax=ax)
#     PF.Simple3DScatter(cores[0,0,gen,i],cores[0,1,gen,i],cores[0,2,gen,i], color=colours[i], marker = 's', fig=fig,ax=ax)

# PF.Simple3DScatter(0,0,0, marker = 'o', color = 'r', fig=fig, ax = ax)

# fig,ax = PF.MakeFig(dimensions='3D')
# gen = 49
# for i in range(0,n_isl):
#     PF.Simple3DScatter(sats_x[0,:,gen,i],sats_y[0,:,gen,i],sats_z[0,:,gen,i], color=colours[i], marker = markers[i], fig=fig,ax=ax)
#     PF.Simple3DScatter(cores[0,0,gen,i],cores[0,1,gen,i],cores[0,2,gen,i], color=colours[i], marker = 's', fig=fig,ax=ax)
    
# PF.Simple3DScatter(0,0,0, marker = 'o', color = 'r', fig=fig, ax = ax)

# #%% 

# file_history = f"propagationHistory_{algo}_sd{sd}_sats{20}_nisl48_npop{n_pop}_int3_d365_best.dat"
# file_core = f"corePosition_{algo}_sd{sd}_sats{5}_nisl48_npop{n_pop}_int3_d365_best.dat"
# # Compute initial satellite positions relative to core
# t, x, y, z, vx, vy, vz,ddd  = DI.ImportPropagationHistory(datafolder + file_history,0, False)
# sat = np.concatenate((x,y,z), axis=1)

# #%% 
# BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
# PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
#                    baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
# #%% 
# N_baseline_too_large = np.sum(BL_m > 100e3)
# N_baseline_too_small = np.sum(BL_m < 500)
# print( f"  this optima has {N_baseline_too_large} instances with too large baselines, \n \
#       and {N_baseline_too_small} instances with too small baselines")