# -*- coding: utf-8 -*-
"""
Animated figure generation for presentation
Created on Mon Oct  5 14:14:18 2020

@author: USER
"""

from FunctionContainer import *

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

#%% Set proper folder for output

figure_folder = LL.folder_dict["Pres_media"]
LL.CheckDir(figure_folder)

#%% Select data to import
cwd = getcwd()
datafolder = cwd + '/Data/Optimization/champions_propagated/'

set_name = '35sat_champ'


yr = int(365.25 * 24 / 4)
#%% 
PF.CloseAll()

# Avoid unnecessarily running over this multiple times
if 'doOnce' not in globals():
    print("run setup!")
    file_history = f"propagationHistory_{set_name}.dat"
    file_core = f"corePosition_{set_name}.dat"
    # compute initial satellite positions relative to core
    t, x, y, z, vx, vy, vz,ddd  = DI.ImportPropagationHistory(datafolder + file_history,0, False)
    sat = np.concatenate((x,y,z), axis=1)
    
    BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
    moon = DI.ImportMoonHistory(datafolder + f"propagationHistory_{set_name}_moon.dat")
    L4 = DI.computeL4Location(moon)
    
    core = DI.ImportSingleColumnFile(datafolder + file_core)
    core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)
    #barycentric satellite motion
    co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
    # barycentric core motion
    core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)
    
    conf_x, conf_y, conf_z = co_x-core_x, co_y-core_y, co_z-core_z
    
    t = (t-t[0])/(24*3600)
    doOnce = True
    
    
oneyrspd = 2
fiveyrspd = 6
print("lets-a-go")

#%%
logic = np.ones(BL_m.shape)
logic[BL_m > 100e3] = 0
logic[BL_m < 500] = 0
logic[BLr_m > 1] = 0

retention = np.sum (logic,axis = 1) /  BL_m.shape[1] * 100


pp.figure(figsize = (15,5))
pp.plot(t/365.25, retention)
pp.xlabel('Time [years]')
pp.ylabel('Viable baseline percentage (%)')
pp.title('Baseline retention over time')
pp.xlim([0,5])
pp.tight_layout()
pp.savefig(figure_folder + f'baseline_retention_{set_name}')

#%% 1yr combined orbit figure

# combinedOrbitfig = combinedOrbitFigure()

# combinedOrbitfig.updateData(core_x[:yr],core_y[:yr],core_z[:yr])

# AnimateDataRangePlot(combinedOrbitfig, f'1yrco_{set_name}', figure_folder + '/1yrCO/',skip=oneyrspd)

# Orbitfig = OrbitFigure()

# Orbitfig.updateData(core_x[:yr],core_y[:yr],core_z[:yr])

# AnimateDataRangePlot(Orbitfig, f'1yrco_xy_{set_name}', figure_folder + '/1yrCO/',skip=oneyrspd)

# #%% 5 yr swarm orbit
# del combinedOrbitfig

# combinedOrbitfig = combinedOrbitFigure()

# combinedOrbitfig.updateData(core_x,core_y,core_z,1)

# AnimateDataRangePlot(combinedOrbitfig, '5yrco', figure_folder + '/5yrCO/',skip=fiveyrspd)

#%% Motion relative to core 1yr
# relToCore = relativeToCoreMotion()
# relToCore.updateData(conf_x[:yr],conf_y[:yr],conf_z[:yr])
# relToCore.pl = 180
# relToCore.updateLims([-50,50],[-50,50],[-50,50])
# AnimateDataRangePlot(relToCore, f'1yrcore_{set_name}', figure_folder + '/1yrCORE/',skip=oneyrspd)

# #%% Motion relative to core 5yr
# relToCore = relativeToCoreMotion()
# relToCore.updateData(conf_x,conf_y,conf_z)
# relToCore.pl = 180
# # relToCore.updateLims([-50,50],[-50,50],[-50,50])
# AnimateDataRangePlot(relToCore, '5yrcore', figure_folder + '/5yrCORE/',skip=fiveyrspd)

#%% Harmonic decay figure 1yr

# harmon = harmonicDecayPlot()
# harmon.updateData(conf_x[:yr],conf_y[:yr], conf_z[:yr], t[:yr])
# AnimateDataRangePlot(harmon, f'1yrharmonics_{set_name}', figure_folder + '/1yrHARM/', skip=oneyrspd)

#%% Harmonic decay figure 5yr

# harmon = harmonicDecayPlot()
# harmon.xlim = [0,5*365]
# harmon.updateData(conf_x,conf_y, conf_z, t)

# AnimateDataRangePlot(harmon, '5yrharmonics', figure_folder + '/5yrHARM/', skip=fiveyrspd)

# #%% Baseline history figure 1yr

# uvw3pnl = uvwBaselinePlot3Panel()
# uvw3pnl.updateData(BL_x[:yr], BL_y[:yr], BL_z[:yr])
# AnimateDataRangePlot(uvw3pnl, '1yruvw3panel', figure_folder + '/1yruvw3panel/', skip=4*oneyrspd)

# #%% Baseline history figure 5yr

# uvw3pnl = uvwBaselinePlot3Panel()
# uvw3pnl.updateData(BL_x, BL_y, BL_z)
# AnimateDataRangePlot(uvw3pnl, '5yruvw3panel', figure_folder + '/5yruvw3panel/', skip=fiveyrspd)

# #%% 3d uvw baseline figure 1yr

# uvw3d = uvwBaselinePlot3D()
# uvw3d.updateData(BL_x[:yr], BL_y[:yr], BL_z[:yr])
# AnimateDataRangePlot(uvw3d, '1yruvw3d', figure_folder + '/1yruvw3d/', skip=4*oneyrspd)

#%% Baseline profile 1yr


# blprof = baselineProfile()
# blprof.updateData(BL_m[:yr], BLr_m[:yr], t[:yr])
# AnimateDataRangePlot(blprof, '1yrblprofile', figure_folder + '/1yrblprofile/', skip=oneyrspd)