# -*- coding: utf-8 -*-
"""
Analyse the results of the coarse optimisation analysis
Created on Sun May 31 14:29:34 2020

@author: Jurriaan
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

figure_folder = LL.folder_dict_funcs["OPT"] + "/Coarse/"
figfolder = LL.folder_dict["NumericalSim"]

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization/Coarse/'

#%% Set up data import loop

PF.ClearPlots()

algorithms = ["de1220"]
algorithms_full = ["Differential Evolution"]

n_sats = [5,10]
n_isl = [1, 4]
n_lcl = [1,5]
n_pop = [128,32]
sd = 42
n_gen = 50

comp_times = [12548,3772,25043,7393, 13330, 35415]
comp_times_compensated = np.empty_like(comp_times)
comp_legends = ['']*6
comp_load, comp_score = np.zeros(6),np.zeros(6)
topgen = np.empty_like(comp_times)

found_null = [False]*6
Npops, Nislands = np.zeros(6),np.zeros(6)
#keep this to easily exchange algorithms using an index
algo = algorithms[0]
algo_name = algorithms_full[0]

# vectors to store the initial configuration coordinates rel to swarm core:
conf_x, conf_y, conf_z = (np.NaN*np.zeros((6,10)) for i in range(0,3))

init = False

plot_colours = [ 'r', 'g', 'b']

#preload moon position data
moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")

c = 0
for N_sats in n_sats:
    for N_isl in n_isl:
        
        N_pop = n_pop[np.where(np.array(n_isl) == N_isl)[0][0]] #not clean but works
                
        file_history = f"propagationHistory_{algo}_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
        file_core = f"corePosition_{algo}_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int1_best.dat"
        for gen in range(1,n_gen+1):
            
            
            comp_load[c] = n_gen*128*1
            
            # also need to loop over the different islands
            for i in range(0,N_isl):
            

                file_fitness = f"fitness_{algo}_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{1}_g{gen}_i{i}.dat" 
                
                fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
                
                if init == False: 
                    fit = np.zeros((len(fitness),n_gen,6,np.max(n_isl)))
                    init=True
                    
                    
                fit[:N_pop,gen-1,c,i] = fitness
                
        comp_legends[c] = f"{algo},{N_isl} islands,{N_pop} population, {1} internal" 
        plot_title = f"{N_sats} satellites, {N_isl} islands, 1 internal"
        file_name = f"{algo}_sd{sd}_sats{N_sats}_isl{N_isl}_pop{N_pop}_int1" 
        fit_mean = np.mean(fit[:N_pop,:,c,:i+1], axis=0)
        fit_max = np.max(fit[:N_pop,:,c,:i+1],axis=0)
        fit_min = np.min(fit[:N_pop,:,c,:i+1],axis=0)
        
        topgen[c] = np.where(np.min(fit_min,axis=1) == 0)[0][0] if (len(np.where(np.min(fit_min,axis=1) == 0)[0]) > 0) else n_gen
        comp_times_compensated[c] = comp_times[c]*(topgen[c]/n_gen)
        
        comp_score[c] = np.sum(np.min(fit_min,axis=1))
        
        Npops[c] = N_pop
        Nislands[c] = N_isl
        
        gens = np.arange(0,n_gen)
        
        fig,ax = PF.MakeFig(figsize=(5,10))
        
        for j in range(0,fit_mean.shape[1]):
            pp.plot(gens,fit_mean[:,j]/1e3,c=plot_colours[0])
            pp.plot(gens,fit_min[:,j]/1e3,c=plot_colours[1])
            pp.plot(gens,fit_max[:,j]/1e3,c=plot_colours[2])

        PF.DressFig(fig, ax, title= plot_title,
                    savefig = True, folder = figure_folder,
                    name = file_name,
                    xlabel = "Generation", ylabel = "Cost [thousands]" ,
                    ylim = [0,np.max(fit_max/1e3)],
                    legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
                    logScale = False, logScaleX = False)
        
        # Compute initial satellite positions relative to core
        t, x, y, z, vx, vy, vz,ddd = DI.ImportPropagationHistory(datafolder + file_history,0, False)
        sat = np.concatenate((x,y,z), axis=1)
    
        core = DI.ImportSingleColumnFile(datafolder + file_core)
        core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)
        co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
        core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)
        
        
        conf_x[c,:N_sats], conf_y[c,:N_sats], conf_z[c,:N_sats] = co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0]
        
        found_null[c] = True if np.any(np.min(fit_min,axis=1) == 0) else False
        
        #%% 
        BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
        PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
                           baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
        PF.saveFig(figfolder ,f"Baselines_{N_sats}_" + comp_legends[c])
        #%% 
        N_baseline_too_large = np.sum(BL_m > 100e3)
        N_baseline_too_small = np.sum(BL_m < 500)
        print( f" {plot_title}: this optima has {N_baseline_too_large} instances with too large baselines, \n \
              and {N_baseline_too_small} instances with too small baselines")
        
        c+=1
        
for N_sats in n_sats:
    comp_load[c] = n_gen*128*5
    
    file_history = f"propagationHistory_{algo}_sd{sd}_sats{N_sats}_nisl4_npop{N_pop}_int5_best.dat"
    file_core = f"corePosition_{algo}_sd{sd}_sats{N_sats}_nisl4_npop{N_pop}_int5_best.dat"
    
    for gen in range(1,n_gen+1):
           
       # also need to loop over the different islands
       for i in range(0,N_isl):
       
           N_pop = n_pop[np.where(np.array(n_isl) == N_isl)[0][0]] #not clean but works
           
           file_fitness = f"fitness_{algo}_sd{sd}_sats{N_sats}_nisl4_npop32_int{5}_g{gen}_i{i}.dat" 

           fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
           
           if init == False: 
               fit = np.zeros((len(fitness),n_gen,6,np.max(n_isl)))
               init=True
               
               
           fit[:N_pop,gen-1,c,i] = fitness
                
    comp_legends[c] = f"{algo},{N_isl} islands,{N_pop} population, {5} internal"    
    plot_title = f"{N_sats} satellites, 4 islands, 5 internal"
    file_name = f"{algo}_sd{sd}_sats{N_sats}_isl{N_isl}_pop{N_pop}_int5" 
    fit_mean = np.mean(fit[:N_pop,:,c,:i+1], axis=0)
    fit_max = np.max(fit[:N_pop,:,c,:i+1],axis=0)
    fit_min = np.min(fit[:N_pop,:,c,:i+1],axis=0)
    
    topgen[c] = np.where(np.min(fit_min,axis=1) == 0)[0][0]
    comp_times_compensated[c] = comp_times[c]*(topgen[c]/n_gen)
    
    comp_score[c] = np.sum(np.min(fit_min,axis=1))
    
    Npops[c] = N_pop
    Nislands[c] = N_isl
    
    gens = np.arange(0,n_gen)
    
    fig,ax = PF.MakeFig(figsize=(5,10))
    
    for j in range(0,fit_mean.shape[1]):
        pp.plot(gens,fit_mean[:,j]/1e3,c=plot_colours[0])
        pp.plot(gens,fit_min[:,j]/1e3,c=plot_colours[1])
        pp.plot(gens,fit_max[:,j]/1e3,c=plot_colours[2])

    PF.DressFig(fig, ax, title= plot_title,
                savefig = True, folder = figure_folder,
                name = file_name,
                xlabel = "Generation", ylabel = "Cost [thousands]" ,
                ylim = [0,np.max(fit_max)/1e3],
                legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
                logScale = False, logScaleX = False)
    
    # Compute initial satellite positions relative to core
    t, x, y, z, vx, vy, vz,ddd  = DI.ImportPropagationHistory(datafolder + file_history,0, False)
    sat = np.concatenate((x,y,z), axis=1)

    core = DI.ImportSingleColumnFile(datafolder + file_core)
    core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)
    co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
    core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)
    
    
    conf_x[c,:N_sats], conf_y[c,:N_sats], conf_z[c,:N_sats] = co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0]

    found_null[c] = True if np.any(np.min(fit_min,axis=1) == 0) else False
    
    
    #%% 
    BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
    PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
                       baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
    PF.saveFig(figfolder ,f"Baselines_{N_sats}_" + comp_legends[c])
    #%% 
    N_baseline_too_large = np.sum(BL_m > 100e3)
    N_baseline_too_small = np.sum(BL_m < 500)
    print( f" {plot_title}: this optima has {N_baseline_too_large} instances with too large baselines, \n \
          and {N_baseline_too_small} instances with too small baselines")
    
    c+=1


        
#%% Evaluate computational load plot 

markers = ['.', 'd', '.', 'd'] + ['^']*2
avg_times = comp_times/comp_load

fig,ax = PF.MakeFig(figsize=(10,10))
for i in range(0,len(comp_times)):
    pp.scatter(comp_score[i], avg_times[i], marker = markers[i])

PF.DressFig(fig,ax, 
            title = 'Average time per evaluation vs. optimisation performance for 50 generations',
            xlabel = 'Sum of champion costs over 50 generations', 
            ylabel = 'Average evaluation time', legendlist = comp_legends)    

    
fig,ax = PF.MakeFig(figsize=(10,10))
for i in range(0,len(comp_times)):
    pp.scatter(comp_times[i],comp_load[i], marker = markers[i])

PF.DressFig(fig,ax, 
            title = 'Function evaluations vs. time required for 50 generations',
            xlabel = 'Computation time [s]', 
            ylabel = 'Function evaluations', legendlist = comp_legends, logScale = True)    


#%% Compensated time to 0

fig,ax = PF.MakeFig(figsize=(10,10))
for i in range(0,len(comp_times)):
    pp.scatter(comp_score[i],comp_times_compensated[i], marker = markers[i])

PF.DressFig(fig,ax, 
            title = 'Compensated total evaluation time vs. average time per evaluation',
            xlabel = 'Sum of champion costs over 50 generations', 
            ylabel = 'Evaluation time till optima', legendlist = comp_legends) 

#%% Comparison of found optima

fig,ax = PF.MakeFig(figsize=(10,10),dimensions = '3D')
for i in range(0,6):
    if found_null[i]: ax.scatter(conf_x[i,:]/1000,conf_y[i,:]/1000,conf_z[i,:]/1000, marker = markers[i], s = 50)


PF.DressFig(fig,ax, 
            title = 'Compensated total evaluation time vs. average time per evaluation',
            xlabel = 'x (rel. to core) [km]', 
            ylabel = 'y (rel. to core) [km]',
            zlabel = 'z (rel. to core) [km]',
            legendlist = comp_legends) 

PF.Animate3DPlotToGIF(fig, ax, "3DScatterOptimal", figure_folder)


#PF.Animate3DPlotToGIF(fig, ax, "Initial swarm configuration", figfolder)

#%% Check out results for seed 72:
sd = 72
init = False
for gen in range(1,n_gen+1):
            
    
    # also need to loop over the different islands
    for i in range(0,N_isl):
    

        file_fitness = f"fitness_{algo}_sd{sd}_sats{5}_nisl{N_isl}_npop{N_pop}_int{1}_g{gen}_i{i}.dat" 
        
        fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
        
        if init == False: 
            fit72 = np.zeros((len(fitness),n_gen,1,N_isl))
            init=True
            
            
        fit72[:,gen-1,0,i] = fitness

plot_title = f"Optimization progress: {algo}, seed:{sd},{5} satellites, {N_isl} islands, {N_pop} population, 1 internal"
file_name = f"Optimization progress_{algo}_sd{sd}_sats{5}_isl{N_isl}_pop{N_pop}_int1" 
fit_mean = np.mean(fit72[:N_pop,:,0,:i+1],axis=0)
fit_max = np.max(fit72[:N_pop,:,0,:i+1], axis=0)
fit_min = np.min(fit72[:N_pop,:,0,:i+1], axis=0)

topgen72 = np.where(np.min(fit_min,axis=1) == 0)[0][0] if (len(np.where(np.min(fit_min,axis=1) == 0)[0]) > 0) else n_gen
comp_times_compensated72 = 2680*(topgen72/n_gen)

comp_score72 = np.sum(np.min(fit_min,axis=1))

gens = np.arange(0,n_gen)

fig,ax = PF.MakeFig(figsize=(10,10))

for j in range(0,fit_mean.shape[1]):
    pp.plot(gens,fit_mean[:,j],c=plot_colours[0])
    pp.plot(gens,fit_min[:,j],c=plot_colours[1])
    pp.plot(gens,fit_max[:,j],c=plot_colours[2])

PF.DressFig(fig, ax, title= plot_title,
            savefig = True, folder = figure_folder,
            name = file_name,
            xlabel = "Generation", ylabel = "Cost" ,
            ylim = [0,np.max(fit_max)],
            legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
            logScale = False, logScaleX = False)



#%% Generate 3d plot comparing changes in rd. seed

fig,ax = PF.MakeFig(figsize=(10,10),dimensions = '3D')
ax.scatter(conf_x[1,:]/1000,conf_y[1,:]/1000,conf_z[1,:]/1000, marker = markers[i], s = 50)

file_history = f"propagationHistory_{algo}_sd{sd}_sats{5}_nisl4_npop{N_pop}_int1_best.dat"
file_core = f"corePosition_{algo}_sd{sd}_sats{5}_nisl4_npop{N_pop}_int1_best.dat"
# Compute initial satellite positions relative to core
t, x, y, z, vx, vy, vz,ddd  = DI.ImportPropagationHistory(datafolder + file_history,0, False)
sat = np.concatenate((x,y,z), axis=1)

core = DI.ImportSingleColumnFile(datafolder + file_core)
core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)

conf_x72, conf_y72, conf_z72 = co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0]

ax.scatter(conf_x72/1000,conf_y72/1000,conf_z72/1000, marker = markers[i], s = 50)


PF.DressFig(fig,ax, 
            title = 'Initial swarm positions relative to core, dependent on seed',
            xlabel = 'x (rel. to core) [km]', 
            ylabel = 'y (rel. to core) [km]',
            zlabel = 'z (rel. to core) [km]',
            legendlist = [ '42', '72']) 


#%% Normalized cpu time per core
comp_times2 = np.append(comp_times, [42794])
comps_normalized = comp_times2*np.array([1,4,1,4,4/5,4/5,1])
#%% generate plot

sats = [5,5,10,10,5,10,15]

pp.figure()
for j in range(0,len(sats)):
    pp.scatter(sats[j],comps_normalized[j])

