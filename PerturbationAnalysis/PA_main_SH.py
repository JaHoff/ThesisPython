# -*- coding: utf-8 -*-
"""
Perturbation analysis main file.
Created on Sat Mar  7 17:20:16 2020

@author: mrfan
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
chdir('PerturbationAnalysis/')
#%% SETTINGS
# Draw 3d plots of motion relative to L4 in barycentric frame, needs large time scales to show effects
barycentricL4MotionComparison = False;

addendumList = ["_1km", "_10km", "_100km"]

addendum = "_100km_SH"
# Number of cases to be taken into account
numcases = 11
# Moon-model affecting border
moonborder = 16


#%% Data import
chdir("../")
cwd = getcwd()
append = cwd + "/Data/"

print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")

# x,y,z,vx,vy,vz, R, V = (np.ndarray([]) for i in range(0,8))

datafolder = cwd + '/Data/PA_SH/'

moon = np.ndarray([])

# Load all the data into 3-dimensional arrays
madevars = False

legend = ["Reference model"]
for j in range(0,3):
    
    for i in range( 1 , numcases):
        filename = datafolder + "propagationHistory_" +str(j)+ "_" + str(i) + addendum + ".dat"
        
        t, xd, yd, zd, vxd, vyd, vzd,mn = DI.ImportPropagationHistory(filename,0, False)
        
        if madevars == False:
            L = xd.shape[0]
            x,y,z,vx,vy,vz, R,V, E  = (np.zeros((xd.shape[0],xd.shape[1],numcases,3)) 
                                       for i in range(0,9))
            Dx,Dy,Dz,Dvx,Dvy,Dvz,Dr,Dv,DDr  = (np.zeros((xd.shape[0],xd.shape[1],numcases,3)) 
                                           for i in range(0,9))
            RelPos = np.zeros((L,3,numcases,3))
            DRelPos = np.zeros((L,3,numcases,3))
            MDRelPos = np.zeros((L,1,numcases,3))
            MDmoon = np.zeros((L,1,numcases,3))
            madevars = True
    
        x[:,:,i,j] = xd
        y[:,:,i,j] = yd
        z[:,:,i,j] = zd
        vx[:,:,i,j] = vxd
        vy[:,:,i,j] = vyd
        vz[:,:,i,j] = vzd
        
        R[:,:,i,j] = np.sqrt(xd**2 + yd**2 + zd**2)
        V[:,:,i,j] = np.sqrt(vxd**2 + vyd**2 + vzd**2)
        
        RelPos[:,:,i,j] = np.concatenate((xd[:,0].reshape(L,1)-xd[:,1].reshape(L,1) ,
                                          yd[:,0].reshape(L,1) - yd[:,1].reshape(L,1),
                                          zd[:,0].reshape(L,1) - zd[:,1].reshape(L,1)),
                                          axis=1)
        
        G = 6.67430e-11
        me = 5.97237E23
        ms = 5
        E[:,:,i,j] = 0.5*ms*V[:,:,i,j]**2 - G*me*ms/R[:,:,i,j]

        Dx[:,:,i,j] = x[:,:,1,0] - xd
        Dy[:,:,i,j] = y[:,:,1,0] - yd
        Dz[:,:,i,j] = z[:,:,1,0] - zd
        Dvx[:,:,i,j] = vx[:,:,1,0] - vxd
        Dvy[:,:,i,j] = vy[:,:,1,0] - vyd
        Dvz[:,:,i,j] = vz[:,:,1,0] - vzd
        
        Dr[:,:,i,j] = np.sqrt(Dx[:,:,i,j]**2 + Dy[:,:,i,j]**2 + Dz[:,:,i,j]**2)
        DDr[:,:,i,j] = Dr[:,:,i,j] - Dr[:,:,i-1,j]
        Dv[:,:,i,j] = np.sqrt(Dvx[:,:,i,j]**2 + Dvy[:,:,i,j]**2 + Dvz[:,:,i,j]**2)
        
        # DE[:,:,i,j] = E[:,:,1,0] - E[:,:,i,j]
    
        DRelPos[:,:,i,j] = RelPos[:,:,1,0] - RelPos[:,:,i,j]
        MDRelPos[:,0,i,j] = np.sqrt(DRelPos[:,0,i,j]**2 + DRelPos[:,1,i,j]**2 
                                    + DRelPos[:,2,i,j]**2)
          
        if j == 0:
            break
        elif j == 1:
            legend.append( "D/O({},{}) - D/O({},{})".format(2*i,2*i,2*(i-1),2*(i-1)))
        elif j == 2:
            legend.append( "D/O({},{}) - D/O({},{})".format(2*i,2*i,2*(i-1),2*(i-1)))
        
        
        
print("Succesfully loaded all data")
del xd, yd, zd, vxd, vyd, vzd

legend = [w.replace("D/O(0,0)","point mass") for w in legend]

t -= t[0]

#%%  Plot preparation
Dr[Dr == 0] = np.NaN

PF.ClearPlots()

styles = [None]*16

#%% Do the plotting

d = LL.folder_dict["NumericalSim"]
#d = LL.folder_dict_funcs["PA"]
markerspacing = 500

PF.Matrix2DPlot(t,DDr[:,:,1:,1], title='Changes in position from Earth SH Degree and orders',
                xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1,ylim=[1e-10,1e4],
                savefig = True, figfolder = d, name='Changes position Earth SH' + addendum,
                legendlist = legend[1:], stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(8,10))

PF.Matrix2DPlot(t,DDr[:,:,1:,2], title='Changes in position from Lunar SH Degree and orders',
                xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1,ylim=[1e-10,1e4],
                savefig = True, figfolder = d, name='Changes position Lunar SH' + addendum,
                legendlist = legend[11:], stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(8,10))
# PF.Matrix2DPlot(t,np.abs(DE[:,:,1:moonborder]), title='(Absolute)Changes in satellite energy from perturbations',
#                 xlabel = 't[days]',ylabel='Total energy [kJ]',xmod=24*3600,
#                 savefig = True, figfolder = d, name='Changes System Energy' + addendum,
#                 legendlist = legend[1:moonborder], stylelist = styles,markN = markerspacing,
#                 logScale = True, figsize=(14,8))

# PF.Matrix2DPlot(t,MDRelPos[:,0,1:moonborder], title='Magnitude change in baseline over time from perturbations',
#                 xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1, ylim=[1e-6, 1e6],
#                 savefig = True, figfolder = d, name='Changes Relative Baseline Satellite Modelling' + addendum,
#                 legendlist = legend[1:moonborder], stylelist = styles[1:], markN = markerspacing,
#                 logScale = True, figsize=(14,8))

