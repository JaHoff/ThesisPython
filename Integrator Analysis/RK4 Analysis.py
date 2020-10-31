# -*- coding: utf-8 -*-
"""
Analysis of integration setting results
Created on Wed Apr  1 14:30:53 2020

@author: Jurriez
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
chdir('IntegratorAnalysis/')

#%% SETTINGS
# Draw 3d plots of motion relative to L4 in barycentric frame, needs large time scales to show effects
barycentricL4MotionComparison = False;

addendum = "_5"
# Number of cases to be taken into account
numcases = 5


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

list_t, list_x, list_y, list_z = ([] for i in range(0,4))
list_Bx, list_By,list_Bz = ([] for i in range(0,3))
for i in range( 0 , numcases):
    filename = datafolder + "propagationHistory_" + str(i) + addendum + ".dat"
    
    t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,1, True)
    list_t.append(t)
    list_x.append(xd)
    list_y.append(yd)
    list_z.append(zd)
    
    if madevars == False:
        Bx, By, Bz = (np.NaN*np.zeros((xd.shape[0],1,numcases))  for i in range(0,3))
        Dx,Dy,Dz,DBx, DBy, DBz, DR, DBR  = (np.NaN*np.zeros((xd.shape[0],2,numcases)) 
                                       for i in range(0,8))
        madevars = True

    list_Bx.append( xd[:,0] - xd[:,1])
    list_By.append( yd[:,0] - yd[:,1])
    list_Bz.append( zd[:,0] - zd[:,1])
    if len(list_t[i]) > len(list_t[0]):
        mask = np.in1d(list_t[i],list_t[0])
        Dx[:,:,i] = list_x[0] - list_x[i][mask,:]
        Dy[:,:,i] = list_y[0] - list_y[i][mask,:]
        Dz[:,:,i] = list_z[0] - list_z[i][mask,:]
        DBx[:,0,i] = list_Bx[0] - list_Bx[i][mask]
        DBy[:,0,i] = list_By[0] - list_By[i][mask]
        DBz[:,0,i] = list_Bz[0] - list_Bz[i][mask]
    else: 
        mask = np.in1d(list_t[0],list_t[i])
        Dx[mask,:,i] = list_x[0][mask,:] - list_x[i]
        Dy[mask,:,i] = list_y[0][mask,:] - list_y[i]
        Dz[mask,:,i] = list_z[0][mask,:] - list_z[i]
        DBx[mask,0,i] = list_Bx[0][mask] - list_Bx[i]
        DBy[mask,0,i] = list_By[0][mask] - list_By[i]
        DBz[mask,0,i] = list_Bz[0][mask] - list_Bz[i]
        
    DR[:,:,i] = np.sqrt(Dx[:,:,i]**2+Dy[:,:,i]**2+Dz[:,:,i]**2)
    DBR[:,0,i] = np.sqrt(DBx[:,0,i]**2+DBy[:,0,i]**2+DBz[:,0,i]**2)
t -= t[0]

#%%  Plot preparation
#Dr[Dr == 0] = np.NaN

PF.clearPlots()

base = float(addendum[1:])
f1 = np.round(base/4,1)
f2 = np.round(base/2,1)
f3 = np.round(base*2,1)
f4 = np.round(base*4,1)

legend = [str(base) + " minute time steps",
          str(f1) +" minute time steps",
          str(f2) + " minute time steps",
          str(f3) +" minute time steps",
          str(f4) +" minute time steps"
    ]

styles = [None, None, None, None, None]

#%% Do the plotting

# Get proper output folder directory to save figures in
d = LL.folder_dict["NumericalSim"]
d = LL.folder_dict_funcs["RK4"]
LL.CheckDir(d)

max1 = np.max(DR[np.isfinite(DR)])
max2 = np.max(DBR[np.isfinite(DBR)])

order = np.max([np.ceil(np.log10(max1)),np.ceil(np.log10(max2))])

ylim = [1e-6, 10**order]

fig, ax = PF.Matrix2DPlot(list_t[0],DR[:,:,1:], title="Changes in satellite positions with respect to "+ addendum[1:]+ " minute step size",
                    xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1,ylim = ylim,
                    savefig = True, figfolder = d, name="Changes Satellite Position" + addendum[1:] + "m",
                    legendlist = legend[1:], stylelist = styles,
                    logScale = True, figsize=(14,8))

fig, ax = PF.Matrix2DPlot(list_t[0],DBR[:,0,1:], title="Changes in satellite baselines with respect to "+ addendum[1:]+ " minute step size",
                    xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1,ylim = ylim,
                    savefig = True, figfolder = d, name="Changes Satellite Baseline" + addendum[1:] + "m",
                    legendlist = legend[1:], stylelist = styles,
                    logScale = True, figsize=(14,8))
