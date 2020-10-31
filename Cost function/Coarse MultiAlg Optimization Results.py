# -*- coding: utf-8 -*-
"""
Created on Sun May 31 14:29:34 2020

@author: mrfan
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

figure_folder = LL.folder_dict["NumericalSim"]

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization/Coarse_Multi_2/'

#%% Set up data import loop

PF.ClearPlots()

algorithms = ["de1220"]
algorithms_full = ["Differential Evolution"]

N_sats = 15
N_isl = 8
n_lcl = 1
N_pop = 32
sd = 72
n_gen = 50

# comp_times = [12548,3772,25043,7393, 13330, 35415]
# comp_times_compensated = np.empty_like(comp_times)
# comp_legends = ['']*6
# comp_load, comp_score = np.zeros(6),np.zeros(6)

# found_null = [False]*6
# Npops, Nislands = np.zeros(6),np.zeros(6)
#keep this to easily exchange algorithms using an index
algo = "de1220"
algo_name = algorithms_full[0]

# vectors to store the initial configuration coordinates rel to swarm core:
conf_x, conf_y, conf_z = (np.NaN*np.zeros((6,10)) for i in range(0,3))

init = False

plot_colours = [ 'r', 'g', 'b']

#preload moon position data
moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")


file_history_single_notop = f"propagationHistory_singlemethod_notop_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
file_core_single_notop = f"corePosition_singlemethod_notop_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
        
file_history_single = f"propagationHistory_singlemethod_top_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
file_core_single = f"corePosition_singlemethod_top_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
file_history_multi = f"propagationHistory_multimethod_top_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
file_core_multi = f"corePosition_multimethod_top_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
for gen in range(1,n_gen+1):    
    # also need to loop over the different islands
    for i in range(0,N_isl):
        file_fitness = f"fitness_singlemethod_top_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{1}_g{gen}_i{i}.dat" 
        fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
        
        
        if init == False: 
            fit_single,fit_multi,fit_notop = np.zeros((len(fitness),n_gen,N_isl)),np.zeros((len(fitness),n_gen,N_isl)), \
                                                np.zeros((len(fitness),n_gen,N_isl))
            init=True
        
        
        fit_single[:N_pop,gen-1,i] = fitness/1e3
        
        file_fitness = f"fitness_singlemethod_notop_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{1}_g{gen}_i{i}.dat" 
        fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
        fit_notop[:N_pop,gen-1,i] = fitness/1e3
        
        file_fitness = f"fitness_multimethod_top_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{1}_g{gen}_i{i}.dat" 
        fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
        fit_multi[:N_pop,gen-1,i] = fitness/1e3

        
plot_title = f"8 {algo} islands, no topology"
file_name = f"8_de1220_no_topology" 
fit_mean = np.mean(fit_notop[:N_pop,:,:i+1], axis=0)
fit_max = np.max(fit_notop[:N_pop,:,:i+1],axis=0)
fit_min = np.min(fit_notop[:N_pop,:,:i+1],axis=0)

topgen = np.where(np.min(fit_min,axis=1) == 0)[0][0] if (len(np.where(np.min(fit_min,axis=1) == 0)[0]) > 0) else n_gen

gens = np.arange(0,n_gen)

fig,ax = PF.MakeFig(figsize=(5,10))

for j in range(0,fit_mean.shape[1]):
    pp.plot(gens,fit_mean[:,j],c=plot_colours[0])
    pp.plot(gens,fit_min[:,j],c=plot_colours[1])
    pp.plot(gens,fit_max[:,j],c=plot_colours[2])

PF.DressFig(fig, ax, title= plot_title,
            savefig = True, folder = figure_folder,
            name = file_name,
            xlabel = "Generation", ylabel = "Cost [thousands]" ,
            ylim = [0,np.max(fit_max)],
            legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
            logScale = False, logScaleX = False)

plot_title = f"8 {algo} islands, full topology"
file_name = f"8_de1220_full_topology" 
fit_mean = np.mean(fit_single[:N_pop,:,:i+1], axis=0)
fit_max = np.max(fit_single[:N_pop,:,:i+1],axis=0)
fit_min = np.min(fit_single[:N_pop,:,:i+1],axis=0)

topgen = np.where(np.min(fit_min,axis=1) == 0)[0][0] if (len(np.where(np.min(fit_min,axis=1) == 0)[0]) > 0) else n_gen

gens = np.arange(0,n_gen)

fig,ax = PF.MakeFig(figsize=(5,10))

for j in range(0,fit_mean.shape[1]):
    pp.plot(gens,fit_mean[:,j],c=plot_colours[0])
    pp.plot(gens,fit_min[:,j],c=plot_colours[1])
    pp.plot(gens,fit_max[:,j],c=plot_colours[2])

PF.DressFig(fig, ax, title= plot_title,
            savefig = True, folder = figure_folder,
            name = file_name,
            xlabel = "Generation", ylabel = "Cost [thousands]" ,
            ylim = [0,np.max(fit_max)],
            legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
            logScale = False, logScaleX = False)

plot_title = f"8 de1220 islands, 1 pso, full topology"
file_name = f"8_multi_full_topology" 
fit_mean = np.mean(fit_multi[:N_pop,:,:i+1], axis=0)
fit_max = np.max(fit_multi[:N_pop,:,:i+1],axis=0)
fit_min = np.min(fit_multi[:N_pop,:,:i+1],axis=0)

topgen = np.where(np.min(fit_min,axis=1) == 0)[0][0] if (len(np.where(np.min(fit_min,axis=1) == 0)[0]) > 0) else n_gen

gens = np.arange(0,n_gen)

fig,ax = PF.MakeFig(figsize=(5,10))

for j in range(0,fit_mean.shape[1]):
    pp.plot(gens,fit_mean[:,j],c=plot_colours[0])
    pp.plot(gens,fit_min[:,j],c=plot_colours[1])
    pp.plot(gens,fit_max[:,j],c=plot_colours[2])

PF.DressFig(fig, ax, title= plot_title,
            savefig = True, folder = figure_folder,
            name = file_name,
            xlabel = "Generation", ylabel = "Cost [thousands]" ,
            ylim = [0,np.max(fit_max)],
            legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
            logScale = False, logScaleX = False)

# Compute initial satellite positions relative to core
# t, x, y, z, vx, vy, vz,ddd = DI.ImportPropagationHistory(datafolder + file_history_single,0, False)
# sat = np.concatenate((x,y,z), axis=1)

# core = DI.ImportSingleColumnFile(datafolder + file_core)
# core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)
# co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
# core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)


# conf_x[c,:N_sats], conf_y[c,:N_sats], conf_z[c,:N_sats] = co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0]

# found_null[c] = True if np.any(np.min(fit_min,axis=1) == 0) else False

# #%% 
# BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
# PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
#                    baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
# #%% 
# N_baseline_too_large = np.sum(BL_m > 100e3)
# N_baseline_too_small = np.sum(BL_m < 500)
# print( f" {plot_title}: this optima has {N_baseline_too_large} instances with too large baselines, \n \
#       and {N_baseline_too_small} instances with too small baselines")
        
        
