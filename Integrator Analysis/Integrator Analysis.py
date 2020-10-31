# -*- coding: utf-8 -*-
"""
Analysis of the effect of different integrators and settings
Created on Wed Apr  1 14:30:53 2020

@author: Jurriez
"""

from os import getcwd
from os import chdir
import sys 
import numpy as np
import matplotlib.pyplot as pp
from adjustText import adjust_text

sys.path.append('../')
chdir('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF
chdir('Integrator Analysis')

#%% SETTINGS
# Draw 3d plots of motion relative to L4 in barycentric frame, needs large time scales to show effects
barycentricL4MotionComparison = False;

addendum = "_365d"
# Number of cases to be taken into account
numcases = 5

prefixes = ["RK4 - 10s",
            "RK4 - 1min",
            "RK4 - 5min",
            "RK4 - 15min",
            "RKF45",
            "RKF56",
            "RKF78",
            "DOPRI87",
            "ABM",
            "BS"]

affixes = ["1E-6",
           "1E-7",
           "1E-8",
           "1E-9",
           "1E-10",
           "1E-11",
           "1E-12",
           "1E-13",
           "1E-14"]

# Generate the list of labels

legend = []
shortlegend = []
timesteps= ["10s", "1min", "5min", "15min"]
for i in range(0,len(prefixes)):
    
    if i <= 3:
        legend.append(prefixes[i])
        shortlegend.append(timesteps[i])
    else:
        if i > 7:
            for j in affixes[:-2]:
                legend.append(prefixes[i]+" tol:"+j)
                shortlegend.append(j)
        else:
            
            for j in affixes:
                legend.append(prefixes[i]+" tol:"+j)
                shortlegend.append(j)

# shortlegend = [s[-5:] for s in legend]            
#%% Data import
chdir("../")
cwd = getcwd()
append = cwd + "/Data/"

print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")

# x,y,z,vx,vy,vz, R, V = (np.ndarray([]) for i in range(0,8))

datafolder = cwd + '/Data/IntegratorAnalysis/'

madevars = False

N = len(prefixes)
NT = len(legend)


cputimes = DI.ImportPropagationHistory(datafolder + "propagationTime" + addendum + ".dat")[0]
numevals = DI.ImportPropagationHistory(datafolder + "functionEvaluations" + addendum + ".dat",1)[0]

count = 0
for i in range( 0 , N):
    if (i<=3):
        filename = datafolder + "interpolatedHistory_" + str(i) + addendum + ".dat"
        t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,0, False)
        
        if madevars == False:
            L = len(t)
            W = xd.shape[1]
            
            # Make appropiate empty variables
            pos_x, pos_y, pos_z, pos_vx, pos_vy, pos_vz = (
                np.zeros((L,W,NT)) for k in range(0,6))
            
            error_x, error_y, error_z, error_vx, error_vy, error_vz = (
                np.zeros((L,W,NT)) for k in range(0,6))
            
            baseline = np.zeros((L,3,NT))
            error_baseline = np.zeros((L,3,NT))
            error_mag = np.zeros((L,W,NT))
            error_mag_baseline = np.zeros((L,W,NT))
            print("Made variables")
            madevars = True
        
        
        
        pos_x[:,:,i], pos_y[:,:,i], pos_z[:,:,i] = xd,yd,zd
        pos_vx[:,:,i], pos_vy[:,:,i], pos_vz[:,:,i] = vxd,vyd,vzd
        
        baseline[:,0,i] = xd[:,1] - xd[:,0]    
        baseline[:,1,i] = yd[:,1] - yd[:,0]
        baseline[:,2,i] = zd[:,1] - zd[:,0]
        
        if i>0:
            error_x[:,:,i], error_y[:,:,i], error_z[:,:,i] = (xd - pos_x[:,:,0],
                                                              yd - pos_y[:,:,0],
                                                              zd - pos_z[:,:,0])
            error_vx[:,:,i], error_vy[:,:,i], error_vz[:,:,i] = (vxd - pos_vx[:,:,0],
                                                                 vyd - pos_vy[:,:,0],
                                                                 vzd - pos_vz[:,:,0])                                                  
            
            error_mag[:,:,i] = np.sqrt( error_x[:,:,i]**2 + error_y[:,:,i]**2 + 
                                 error_z[:,:,i]**2)
            
            error_baseline[:,:,i]  = baseline[:,:,i] - baseline[:,:,0]
            
            error_mag_baseline[:,0,i] = np.linalg.norm(error_baseline[:,:,i],axis = 1)
            
        count += 1
        
    else:
        if i > 7:
            for j in range(0,len(affixes)-2):
                affix = affixes[j][1]+affixes[j][3:]
                filename = datafolder + "interpolatedHistory_" + str(i) + addendum + "_"+ affix +".dat"
                ij = count
                t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,0, False)
                
                baseline[:,0,ij] = xd[:,1] - xd[:,0]    
                baseline[:,1,ij] = yd[:,1] - yd[:,0]
                baseline[:,2,ij] = zd[:,1] - zd[:,0]
                
                pos_x[:,:,ij], pos_y[:,:,ij], pos_z[:,:,ij] = xd,yd,zd
                pos_vx[:,:,ij], pos_vy[:,:,ij], pos_vz[:,:,ij] = vxd,vyd,vzd
            
    
                error_x[:,:,ij], error_y[:,:,ij], error_z[:,:,ij] = (xd - pos_x[:,:,0],
                                                                     yd - pos_y[:,:,0],
                                                                     zd - pos_z[:,:,0])
                error_vx[:,:,ij], error_vy[:,:,ij], error_vz[:,:,ij] = (vxd - pos_vx[:,:,0],
                                                                        vyd - pos_vy[:,:,0],
                                                                        vzd - pos_vz[:,:,0])                                                  
                
                error_mag[:,:,ij] = np.sqrt( error_x[:,:,ij]**2 + error_y[:,:,ij]**2 + 
                                             error_z[:,:,ij]**2)
                
                error_baseline[:,:,ij]  = baseline[:,:,ij] - baseline[:,:,0]
                
                error_mag_baseline[:,0,ij] = np.linalg.norm(error_baseline[:,:,ij],axis = 1)
                
                count += 1
        else:        
            for j in range(0,len(affixes)):
                affix = affixes[j][1]+affixes[j][3:]
                filename = datafolder + "interpolatedHistory_" + str(i) + addendum + "_"+ affix +".dat"
                ij = count
                t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,0, False)
                
                baseline[:,0,ij] = xd[:,1] - xd[:,0]    
                baseline[:,1,ij] = yd[:,1] - yd[:,0]
                baseline[:,2,ij] = zd[:,1] - zd[:,0]
                
                pos_x[:,:,ij], pos_y[:,:,ij], pos_z[:,:,ij] = xd,yd,zd
                pos_vx[:,:,ij], pos_vy[:,:,ij], pos_vz[:,:,ij] = vxd,vyd,vzd
            
    
                error_x[:,:,ij], error_y[:,:,ij], error_z[:,:,ij] = (xd - pos_x[:,:,0],
                                                                     yd - pos_y[:,:,0],
                                                                     zd - pos_z[:,:,0])
                error_vx[:,:,ij], error_vy[:,:,ij], error_vz[:,:,ij] = (vxd - pos_vx[:,:,0],
                                                                        vyd - pos_vy[:,:,0],
                                                                        vzd - pos_vz[:,:,0])                                                  
                
                error_mag[:,:,ij] = np.sqrt( error_x[:,:,ij]**2 + error_y[:,:,ij]**2 + 
                                             error_z[:,:,ij]**2)
                
                error_baseline[:,:,ij]  = baseline[:,:,ij] - baseline[:,:,0]
                
                error_mag_baseline[:,0,ij] = np.linalg.norm(error_baseline[:,:,ij],axis = 1)
                
                count += 1
                
    

# Correct for initial time to make the time axis relative
t -= t[0]

#%%  Plot preparation

# Get proper output folder directory to save figures in
dp = LL.folder_dict["NumericalSim"]
d = LL.folder_dict_funcs["IA"]
LL.CheckDir(d)

PF.ClearPlots()

# Adjust the affixes for plot labelling
affixes = [" 1E-6",
           " 1E-7",
           " 1E-8",
           " 1E-9",
           "1E-10",
           "1E-11",
           "1E-12"]
#%% Scatter plots

colormap = ['r']*4 + ['b']* 9 + ['g']*9 + ['k']*9 + ['c']*9 + ['m']*7 + ['y']*7
legendentries = ['RK4'] + [None]*3 + ['RKF45'] + [None] * 8 + ['RKF56'] + [None]*8 + \
    ['RKF78'] + [None]*8 +  ['DOPRI87'] + [None]*8 +  ['ABM'] + [None]*6 + ['BS'] + [None]*6    
    
    
    
legend = [w[-5:] for w in legend]
legend = [w.replace(":", " ") for w in legend]

ylims = [1e-4, 1e7]
fig,ax = PF.MakeFig(figsize = (14,8))
for i in range(0,NT):
    PF.Simple2DScatter(cputimes[i], error_mag[-1,0,i],fig = fig, ax = ax ,
                       xmod = 1, ymod = 1, addAnnotation = True, 
                       annotation = legend[i][-5:], color = colormap[i],
                       legendAddition = legendentries[i])
ax.axvline(cputimes[0])
ax.text(cputimes[0]+5, 1, "Reference cpu time")
PF.DressFig(fig, ax, title='Cumulative orbit propagation error versus cpu time',
            xlabel='cpu time [s]', ylabel= 'position error [m]',
            xlim =[4E-2,1E3], ylim = ylims, logScale = True, logScaleX = True)
PF.saveFig(d, "Pos error v cpu time")

fig,ax = PF.MakeFig(figsize = (14,8))
for i in range(0,NT):
    PF.Simple2DScatter(cputimes[i], error_mag_baseline[-1,0,i],fig = fig, ax = ax ,
                       xmod = 1, ymod = 1, addAnnotation = True, 
                       annotation = legend[i][-5:], color = colormap[i],
                       legendAddition = legendentries[i])
ax.axvline(cputimes[0])
ax.text(cputimes[0]+5, 1, "Reference cpu time")
PF.DressFig(fig, ax, title='Cumulative baseline error versus cpu time',
            xlabel='cpu time [s]', ylabel= 'baseline error [m]',
            xlim =[4E-2,1E3], ylim = ylims, logScale = True, logScaleX = True)
PF.saveFig(d, "baseline error v cpu time")

fig,ax = PF.MakeFig(figsize = (14,8))
for i in range(0,NT):
    PF.Simple2DScatter(error_mag[-1,0,i],error_mag_baseline[-1,0,i], fig = fig, ax = ax ,
                       xmod = 1, ymod = 1, addAnnotation = True, 
                       annotation = legend[i][-5:], color = colormap[i],
                       legendAddition = legendentries[i])
ax.axvline(numevals[0])
ax.text(numevals[0]+5, 1, "Reference number of computations")
PF.DressFig(fig, ax, title='Cumulative baseline error over position error',
            xlabel='Position error [m]', ylabel= 'Baseline error [m]',
            xlim =ylims, ylim = ylims, logScale = True, logScaleX = True)
PF.saveFig(d, "Baseline error v pos error")

#%% Better plot for report
textt = []
ScatY = np.max(error_mag_baseline[:,0,:], axis=0)
fig,ax = PF.MakeFig(figsize = (14,8))
for i in range(0,NT):
    a,b,textt = PF.Simple2DScatterAdjustedText(cputimes[i], ScatY[i],fig = fig, ax = ax ,
                       xmod = 1, ymod = 1, addAnnotation = True, 
                       annotation = legend[i], color = colormap[i],
                       legendAddition = legendentries[i],texts = textt)

ax.axvline(cputimes[0])
ax.text(cputimes[0]+5, 1, "Reference cpu time")

pp.axhline(0.1, color = 'r')
ax.text(0.1, 0.11, "Error tolerance")
PF.DressFig(fig, ax, title='Maximum baseline error compared to reference versus cpu time',
            xlabel='cpu time [s]', ylabel= 'baseline error [m]',
            xlim =[4E-2,1E3], ylim = [1e-3,1e4], logScale = True, logScaleX = True)
PF.saveFig(dp, "baseline error v cpu time")

#%% Matrix plots

fig,ax = PF.MakeFig(figsize = (14,8))
PF.Matrix2DPlot(t,error_mag, fig = fig, ax = ax, xmod=24*3600, ymod = 1 ,stylelist = colormap)
PF.DressFig(fig,ax, title= "Total position error over time", xlabel="Propagation time [days]",
            ylabel = "Position error [m]", xlim = [0,365], ylim = ylims,
            logScale = True)

#%% Subfigure plots position
ylims = [1e-6, 5e4]
fig,ax = pp.subplots(1,4, figsize=(12,8))
i = 0
e = 4
PF.Matrix2DPlot(t,error_mag[:,:,i:e], fig = fig, ax = ax[0], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[0:4], stylelist = None)
PF.DressFig(fig,ax[0], title= "RK4 error", xlabel="Propagation time [days]",
            ylabel = "Position error [m]", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e
e = 9

PF.Matrix2DPlot(t,error_mag[:,:,i:i+e], fig = fig, ax = ax[1], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[1], title= "RKF45 error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)
i += e

PF.Matrix2DPlot(t,error_mag[:,:,i:i+e], fig = fig, ax = ax[2], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[2], title= "RKF56 error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e

PF.Matrix2DPlot(t,error_mag[:,:,i:i+e], fig = fig, ax = ax[3], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[3], title= "RKF78 error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e

PF.saveFig(dp, "Position error set 1")
fig.show()

fig,ax = pp.subplots(1,3, figsize=(12,8))

PF.Matrix2DPlot(t,error_mag[:,:,i:i+e], fig = fig, ax = ax[0], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[0], title= "DOPRI87 error", xlabel="Propagation time [days]",
            ylabel = "Position error [m]", xlim = [0,365], ylim = ylims,
            logScale = True)
e = 7
i += e

PF.Matrix2DPlot(t,error_mag[:,:,i:i+e], fig = fig, ax = ax[1], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[1], title= "ABM error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)
i += e

PF.Matrix2DPlot(t,error_mag[:,:,i:i+e], fig = fig, ax = ax[2], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[2], title= "BS error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)

PF.saveFig(dp, "Position error set 2")
fig.show()

#%% Subfigure plots basleine

fig,ax = pp.subplots(1,4, figsize=(16,8))
i = 0
e = 4
PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:e], fig = fig, ax = ax[0], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[0:4], stylelist = None)
PF.DressFig(fig,ax[0], title= "RK4 error", xlabel="Propagation time [days]",
            ylabel = "Baseline error [m]", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e
e = 9

PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:i+e], fig = fig, ax = ax[1], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[1], title= "RKF45 error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)
i += e

PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:i+e], fig = fig, ax = ax[2], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[2], title= "RKF56 error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e

PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:i+e], fig = fig, ax = ax[3], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[3], title= "RKF78 error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e

PF.saveFig(dp, "Baseline error set 1 ")
fig.show()

fig,ax = pp.subplots(1,3, figsize=(12,8))

PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:i+e], fig = fig, ax = ax[0], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[0], title= "DOPRI87 error", xlabel="Propagation time [days]",
            ylabel = "Baseline error [m]", xlim = [0,365], ylim = ylims,
            logScale = True)

i += e
e = 7

PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:i+e], fig = fig, ax = ax[1], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[1], title= "ABM error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)
i += e

PF.Matrix2DPlot(t,error_mag_baseline[:,:,i:i+e], fig = fig, ax = ax[2], xmod=24*3600, ymod = 1 ,
                legendlist = shortlegend[i:i+e], stylelist = None)
PF.DressFig(fig,ax[2], title= "BS error", xlabel="Propagation time [days]",
            ylabel = "", xlim = [0,365], ylim = ylims,
            logScale = True)

PF.saveFig(dp, "Baseline error set 2")
fig.show()
