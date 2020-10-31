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

addendum = "_100km"
# Number of cases to be taken into account
numcases = 25
# Moon-model affecting border
moonborder = 15


#%% Data import
chdir("../")
cwd = getcwd()
append = cwd + "/Data/"

print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")

# x,y,z,vx,vy,vz, R, V = (np.ndarray([]) for i in range(0,8))

datafolder = cwd + '/Data/PerturbationAnalysis/'

moon = np.ndarray([])

# Load all the data into 3-dimensional arrays
madevars = False
for i in range( 0 , numcases):
    filename = datafolder + "propagationHistory_" + str(i) + addendum + ".dat"
    
    t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,1, True)
    
    if madevars == False:
        L = xd.shape[0]
        moon = np.zeros((md.shape[0],md.shape[1],numcases))
        Dmoon = np.zeros((md.shape[0],md.shape[1],numcases))
        x,y,z,vx,vy,vz, R,V, E  = (np.zeros((xd.shape[0],xd.shape[1],numcases)) 
                                   for i in range(0,9))
        Dx,Dy,Dz,Dvx,Dvy,Dvz,Dr,Dv,DE  = (np.zeros((xd.shape[0],xd.shape[1],numcases)) 
                                       for i in range(0,9))
        RelPos = np.zeros((L,3,numcases))
        DRelPos = np.zeros((L,3,numcases))
        MDRelPos = np.zeros((L,1,numcases))
        MDmoon = np.zeros((L,1,numcases))
        madevars = True

    moon[:,:,i] = md
    x[:,:,i] = xd
    y[:,:,i] = yd
    z[:,:,i] = zd
    vx[:,:,i] = vxd
    vy[:,:,i] = vyd
    vz[:,:,i] = vzd
    
    R[:,:,i] = np.sqrt(xd**2 + yd**2 + zd**2)
    V[:,:,i] = np.sqrt(vxd**2 + vyd**2 + vzd**2)
    
    RelPos[:,:,i] = np.concatenate((xd[:,0].reshape(L,1)-xd[:,1].reshape(L,1) ,
                                    yd[:,0].reshape(L,1) - yd[:,1].reshape(L,1),
                                    zd[:,0].reshape(L,1) - zd[:,1].reshape(L,1)),
                                   axis=1)
    
    G = 6.67430e-11
    me = 5.97237E23
    ms = 5
    E[:,:,i] = 0.5*ms*V[:,:,i]**2 - G*me*ms/R[:,:,i]
     
    if i > 0:
        Dmoon[:,:,i] = moon[:,:,0] - moon[:,:,i]
        Dx[:,:,i] = x[:,:,0] - x[:,:,i]
        Dy[:,:,i] = y[:,:,0] - y[:,:,i]
        Dz[:,:,i] = z[:,:,0] - z[:,:,i]
        Dvx[:,:,i] = vx[:,:,0] - vx[:,:,i]
        Dvy[:,:,i] = vy[:,:,0] - vy[:,:,i]
        Dvz[:,:,i] = vz[:,:,0] - vz[:,:,i]
        
        Dr[:,:,i] = np.sqrt(Dx[:,:,i]**2 + Dy[:,:,i]**2 + Dz[:,:,i]**2)
        Dv[:,:,i] = np.sqrt(Dvx[:,:,i]**2 + Dvy[:,:,i]**2 + Dvz[:,:,i]**2)
        
        DE[:,:,i] = E[:,:,0] - E[:,:,i]
        
        DRelPos[:,:,i] = RelPos[:,:,0] - RelPos[:,:,i]
        MDRelPos[:,0,i] = np.sqrt(DRelPos[:,0,i]**2 + DRelPos[:,1,i]**2 + DRelPos[:,2,i]**2)
        
        MDmoon[:,0,i] = np.sqrt(Dmoon[:,0,i]**2 + Dmoon[:,1,i]**2 + Dmoon[:,2,i]**2)
        
print("Succesfully loaded all data")
del xd, yd, zd, vxd, vyd, vzd

t -= t[0]

#%%  Plot preparation
Dr[Dr == 0] = np.NaN

PF.ClearPlots()

legend = ['Default model',
          'Mars point mass removed',
          'Venus point mass removed',
          'Solar point mass removed',
          'Radiation pressure removed',
          'Jupiter point mass removed',
          'Saturn point mass removed',
          'Neptune point mass removed',
          'Uranus point mass removed',
          'Earth as point mass',
          'Earth Harmonics to D/O 10,10',
          'Earth Harmonics to D/O 20,20',
          'Moon harmonics D/O 5,5',
          'Moon harmonics D/O 10,10',
          'Moon harmonics D/O 20,20',
          'Mars point mass removed',
          'Venus point mass removed',
          'Solar point mass removed',
          'Radiation pressure removed',
          'Jupiter point mass removed',
          'Saturn point mass removed',
          'Neptune point mass removed',
          'Uranus point mass removed',
          'Earth point mass',
          'Earth Harmonics D/O 10,10',
          'Earth Harmonics D/O 20,20',
    ]

styles = [None, None, None, '-d', '-d', None, None, None, None, '-o','-o','-o','--', '--', '--']

#%% Do the plotting

d = LL.folder_dict["NumericalSim"]
#d = LL.folder_dict_funcs["PA"]
markerspacing = 500

PF.Matrix2DPlot(t,Dr[:,:,1:moonborder], title='Changes in satellite positions as range over time',
                xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1,
                savefig = True, figfolder = d, name='Changes Satellite Modelling' + addendum,
                legendlist = legend[1:moonborder], stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(14,8))

PF.Matrix2DPlot(t,Dr[:,:,moonborder:], title='Changes in satellite positions from Lunar perturbations as range over time',
                savefig = True, figfolder = d, name='Changes Lunar Modelling' + addendum,
                xlabel = 't[days]',ylabel='change [m]',xmod=24*3600,ymod=1,
                legendlist = legend[moonborder:],stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(14,8))

PF.Matrix2DPlot(t,np.abs(DE[:,:,1:moonborder]), title='(Absolute)Changes in satellite energy from perturbations',
                xlabel = 't[days]',ylabel='Total energy [kJ]',xmod=24*3600,
                savefig = True, figfolder = d, name='Changes System Energy' + addendum,
                legendlist = legend[1:moonborder], stylelist = styles,markN = markerspacing,
                logScale = True, figsize=(14,8))

PF.Matrix2DPlot(t,MDRelPos[:,0,1:moonborder], title='Magnitude change in baseline over time from perturbations',
                xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1, ylim=[1e-6, 1e5],
                savefig = True, figfolder = d, name='Changes Relative Baseline Satellite Modelling' + addendum,
                legendlist = legend[1:moonborder], stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(14,8))

PF.Matrix2DPlot(t,MDRelPos[:,0,moonborder:], title='Magnitude change in baseline over time from perturbations in Lunar modelling',
                xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1, ylim=[1e-6, 1e5],
                savefig = True, figfolder = d, name='Changes Relative Baseline Lunar Modelling' + addendum,
                legendlist = legend[moonborder:], stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(14,8))

PF.Matrix2DPlot(t,MDmoon[:,0,moonborder:], title='Magnitude change in moon position from removal of perturbations',
                xlabel = 't[days]',ylabel='change [m]', xmod=24*3600,ymod = 1, ylim=[1e-6, 5e8],
                savefig = True, figfolder = d, name='Changes Moon position' + addendum,
                legendlist = legend[moonborder:], stylelist = styles[1:], markN = markerspacing,
                logScale = True, figsize=(14,8))

#%% Make scatter plots for the report
scatX = Dr[-1,0,:]
scatY = MDRelPos[-1,0,:]

# Alter the legend list to remove problematic entries
legend = ['Default model',
          'Mars point mass removed',
          'Venus point mass removed',
          'Solar point mass removed',
          'Radiation pressure removed',
          'Jupiter point mass removed',
          'Saturn point mass removed',
          'Neptune point mass removed',
          'Uranus point mass removed',
          'Earth as point mass',
          '',
          '',
          'Moon harmonic models',
          '',
          '',
          'Mars point mass removed',
          'Venus point mass removed',
          'Solar point mass removed',
          '',
          'Jupiter point mass removed',
          'Saturn point mass removed',
          'Neptune point mass removed',
          'Uranus point mass removed',
          'Earth point mass',
          '',
          'Earth Harmonics D/O 20,20',
    ]


texts = []
fig, ax = PF.Simple2DScatter(np.zeros((1,1)) ,np.zeros((1,1)),logScale = True, figsize=(14,8))
for i in range(1,moonborder):
        a,b, texts = PF.Simple2DScatterAdjustedText(scatX[i],scatY[i],title='Changes in relative baseline against absolute position from perturbating forces',
                               xlabel = 'Final change in absolute position (#1) [km]', 
                               ylabel = 'Final change in baseline [km]', xmod = 1.,
                               ymod = 1., xlim = [5E-2,1E10], ylim = [1E-5,1E6] , fig = fig,
                               ax = ax, addAnnotation = True, annotation = legend[i] ,texts = texts)

pp.axhline(0.1, color='r') 
pp.axvline(10, color='r')
pp.text(1e6, 0.2, "10 cm baseline uncertainty limit" , c = 'r'  )
pp.text(1.1e1, 1e3, "10m position uncertainty limit", c='r'  )
adjust_text(texts)
PF.saveFig(d, 'Relative versus absolute' + addendum)

texts = []
fig, ax = PF.Simple2DScatter(np.zeros((1,1)) ,np.zeros((1,1)),logScale = True, figsize=(14,8))
for i in range(moonborder,numcases):
        a,b, texts = PF.Simple2DScatterAdjustedText(scatX[i],scatY[i],title='Changes in relative baseline against absolute position from perturbating forces acting on Moon',
                           xlabel = 'Final change in absolute position (#1) [km]', 
                           ylabel = 'Final change in baseline [km]', xmod = 1.,
                           ymod = 1., xlim = [5E-2,1E10], ylim = [1E-5,1E6] , fig = fig,
                           ax = ax, addAnnotation = True, annotation = legend[i], texts = texts)

pp.axhline(0.1, color='r')
pp.axvline(10, color='r')
pp.text(1e6, 0.2, "10 cm baseline uncertainty limit" , c = 'r'  )
pp.text(1.1e1, 1e3, "10m position uncertainty limit", c='r'  )
adjust_text(texts)
PF.saveFig(d, 'Relative versus absolute Lunar' + addendum)

if barycentricL4MotionComparison:
    for i in range(0,numcases):
        co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon[:,:,i],
                                                              np.concatenate((x[:,:,i],y[:,:,i],
                                                              z[:,:,i]), axis=1),
                                                               exportL4 = True)
    
        N = co_x.shape[0]
        fig2, ax2 = PF.Simple3DPlot(co_x - L4b[:,0].reshape(N,1)*np.ones(co_x.shape), 
                                    co_y - L4b[:,1].reshape(N,1)*np.ones(co_y.shape),
                                    co_z - L4b[:,2].reshape(N,1)*np.ones(co_z.shape))
        PF.Simple3DScatter(0,0,0, fig= fig2, ax = ax2, title=legend[i])
             
#%% Plot the maximum values

Dr[np.isnan(Dr)] = 0
scatX = np.max(Dr[:,0,:],axis = 0)
scatY = np.max(MDRelPos[:,0,:], axis = 0)
Dr[Dr == 0] = np.nan

texts = []
fig, ax = PF.Simple2DScatter(np.zeros((1,1)) ,np.zeros((1,1)),logScale = True, figsize=(14,8))
for i in range(1,moonborder):
        a,b, texts = PF.Simple2DScatterAdjustedText(scatX[i],scatY[i],title='Maximum change in relative baseline against absolute position from perturbating forces',
                               xlabel = 'Final change in absolute position (#1) [km]', 
                               ylabel = 'Final change in baseline [km]', xmod = 1.,
                               ymod = 1., xlim = [5E-2,1E10], ylim = [1E-3,1E6] , fig = fig,
                               ax = ax, addAnnotation = True, annotation = legend[i] ,texts = texts)

pp.axhline(0.1, color='r') 
pp.axvline(10, color='r')
pp.text(1e6, 0.2, "10 cm baseline uncertainty limit" , c = 'r'  )
pp.text(1.1e1, 1e3, "10m position uncertainty limit", c='r'  )
adjust_text(texts)
PF.saveFig(d, 'Maximum relative versus absolute' + addendum)

texts = []
fig, ax = PF.Simple2DScatter(np.zeros((1,1)) ,np.zeros((1,1)),logScale = True, figsize=(14,8))
for i in range(moonborder,numcases):
        a,b, texts = PF.Simple2DScatterAdjustedText(scatX[i],scatY[i],
                            title='Maximum change in relative baseline against absolute position from perturbating forces acting on Moon',
                            xlabel = 'Final change in absolute position (#1) [km]', 
                            ylabel = 'Final change in baseline [km]', xmod = 1.,
                            ymod = 1., xlim = [5E-2,1E10], ylim = [1E-3,1E6] , fig = fig,
                            ax = ax, addAnnotation = True, annotation = legend[i], texts = texts)

pp.axhline(0.1, color='r')
pp.axvline(10, color='r')
pp.text(1e6, 0.2, "10 cm baseline uncertainty limit" , c = 'r'  )
pp.text(1.1e1, 1e3, "10m position uncertainty limit", c='r'  )
adjust_text(texts)
PF.saveFig(d, 'Maximum relative versus absolute Lunar' + addendum)
#%% Small test for linearity of datasets
        
#check position change in x, y and z for 1, 10 and 100 km at time points.