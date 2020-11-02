# -*- coding: utf-8 -*-
"""
Analyse the result of a single coarse optimisation 
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
# import MayaviPlotting as MF
import DataImport as DI
import GeneralFunctions as GF

figure_folder = LL.folder_dict_funcs["OPT"] + "/Coarse/"

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization/Coarse/'

#%% Set up data import loop

PF.ClearPlots()

algorithms = ["de1220"]
algorithms_full = ["Differential Evolution"]

N_sats = 20
N_isl = 5
N_lcl = 5
N_pop = 32
sd = 72
n_gen = 50
algo = "de1220"

# vectors to store the initial configuration coordinates rel to swarm core:
conf_x, conf_y, conf_z = (np.NaN*np.zeros((1,N_sats)) for i in range(0,3))

init = False

plot_colours = [ 'r', 'g', 'b']

#preload moon position data
moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")

c = 0
        
file_history = f"propagationHistory_{algo}_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{N_lcl}_d365_best.dat"
file_core = f"corePosition_{algo}_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{N_lcl}_d365_best.dat"

for gen in range(1,n_gen+1):    
    # also need to loop over the different islands
    for i in range(0,N_isl):
        file_fitness = f"fitness_{algo}_sd{sd}_sats{N_sats}_nisl{N_isl}_npop{N_pop}_int{N_lcl}_d365_g{gen}_i{i}.dat" 
        fitness = DI.ImportSingleColumnFile(datafolder + file_fitness)
        
        if init == False: 
            fit = np.zeros((len(fitness),n_gen,6,N_isl))
            init=True
            
            
        fit[:N_pop,gen-1,c,i] = fitness
        

plot_title = f"Single algorithm: {algo}, seed:{sd},{N_sats} satellites, {N_isl} islands, {N_pop} population, 1 internal"
file_name = f"Single algo_{algo}_sd{sd}_sats{N_sats}_isl{N_isl}_pop{N_pop}_int1" 
fit_mean = np.mean(fit[:N_pop,:,c,:i+1], axis=0)
fit_max = np.max(fit[:N_pop,:,c,:i+1],axis=0)
fit_min = np.min(fit[:N_pop,:,c,:i+1],axis=0)

topgen = np.where(np.min(fit_min,axis=1) == 0)[0][0] if (len(np.where(np.min(fit_min,axis=1) == 0)[0]) > 0) else n_gen

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


#%% Compute initial satellite positions relative to core
t, x, y, z, vx, vy, vz,ddd = DI.ImportPropagationHistory(datafolder + file_history,0, False)
sat = np.concatenate((x,y,z), axis=1)

core = DI.ImportSingleColumnFile(datafolder + file_core)
core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)

l_moon = moon.shape[0]
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat[:l_moon,:], exportL4 = True, fixedL4 = True)
core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos[:l_moon,:], exportL4 = True, fixedL4 = True)


conf_x[:N_sats], conf_y[:N_sats], conf_z[:N_sats] = co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0]


#%% 
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
                    baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
#%% 
N_baseline_too_large = np.sum(BL_m > 100e3)/2
N_baseline_too_small = np.sum(BL_m < 500)/2
vlim = 1
N_baseline_rate_too_large = np.sum(BLr_m > vlim)/2
print( f" {plot_title}: this optima has {N_baseline_too_large} instances with too large baselines, \n \
      and {N_baseline_too_small} instances with too small baselines. \n \
          {N_baseline_rate_too_large} instances where the baseline rate surpassed {vlim} m/s")
      
# fig = MF.MakeFig("movement relative to swarm core",theme="homeworld")
# MF.Plot3D(fig, co_x-core_x, co_y-core_y,
#                           co_z-core_z, lw = 0.1, col = (1,0,0), scale=1e3, opacity = 0.5)

# MF.PlotColorSphere(fig, 3,c=(1,1,0),opacity = 1)
# # MF.AddDefaultAxes(fig,xlabel= 'x km', ylabel= 'y km', zlabel = 'z km')
# MF.PlotAnimatedParticles(fig,co_x-core_x, co_y-core_y,
#                           co_z-core_z, R = 1e3,theme="default", lw=0.2,speed=3)

#%% First swarm config


fig,ax = PF.MakeFig(None,'3D')
PF.Simple3DScatter(0, 0, 0, fig=fig, ax=ax)
PF.Simple3DScatter(co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0],
                   fig=fig, ax=ax,marker = 'd')

PF.Simple3DScatter(L4b[0,0]-core_x[0], L4b[0,1]-core_y[0], L4b[0,2]-core_z[0],
                   fig=fig, ax=ax,marker = '*')
PF.DressFig(fig, ax, 'Intitial swarm configuration relative to swarm core', 
            xlabel = 'km', ylabel = 'km', zlabel = 'km', legendlist = [ "core", "satellites", "L4"])



fig,ax = PF.MakeFig(None,'3D')

PF.Simple3DScatter(L4b[0,0], L4b[0,1], L4b[0,2],
                   fig=fig, ax=ax,marker = '*')

PF.Simple3DPlot(core_x,core_y,core_z, fig=fig, ax = ax)
PF.Simple3DScatter(co_x[0,:], co_y[0,:], co_z[0,:],
                   fig=fig, ax=ax,marker = 'd')
PF.DressFig(fig, ax, 'Trajectory of satellite swarm in barycentric frame', 
            xlabel = 'km', ylabel = 'km', zlabel = 'km', legendlist = [ "Swarm trajectory", "L4","Initial satellite distribution"])
        
#%%
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x[:1,:],y[:1,:],z[:1,:], vx[:1,:], vy[:1,:], vz[:1,:], t[:1])

fig,ax = PF.Simple3DScatter(BL_x.T, BL_y.T, BL_z.T)
PF.Simple3DScatter(0,0,0, marker = 'd',fig = fig, ax = ax)