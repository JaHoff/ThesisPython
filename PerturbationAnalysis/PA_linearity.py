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

import statsmodels.api as sm

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

addendumList = ["_1km", "_10km", "_30km", "_50km", "_70km", "_100km"]

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

for j in range(0, len(addendumList)):
    addendum = addendumList[j]
    
    for i in range( 0 , 5):
        filename = datafolder + "propagationHistory_" + str(i) + addendum + ".dat"
        
        t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,1, True)
        
        if madevars == False:
            L = xd.shape[0]
            moon = np.zeros((md.shape[0],md.shape[1],numcases,len(addendumList)))
            Dmoon = np.zeros((md.shape[0],md.shape[1],numcases,len(addendumList)))
            x,y,z,vx,vy,vz, R,V  = (np.zeros((xd.shape[0],xd.shape[1],numcases,len(addendumList))) 
                                       for i in range(0,8))
            Dx,Dy,Dz,Dvx,Dvy,Dvz,Dr,Dv  = (np.zeros((xd.shape[0],xd.shape[1],numcases,len(addendumList))) 
                                           for i in range(0,8))
            RelPos = np.zeros((L,3,numcases,len(addendumList)))
            DRelPos = np.zeros((L,3,numcases,len(addendumList)))
            MDRelPos = np.zeros((L,1,numcases,len(addendumList)))
            MDmoon = np.zeros((L,1,numcases,len(addendumList)))
            madevars = True
    
        moon[:,:,i,j] = md
        x[:,:,i,j] = xd
        y[:,:,i,j] = yd
        z[:,:,i,j] = zd
        vx[:,:,i,j] = vxd
        vy[:,:,i,j] = vyd
        vz[:,:,i,j] = vzd
        
        R[:,:,i,j] = np.sqrt(xd**2 + yd**2 + zd**2)
        V[:,:,i,j] = np.sqrt(vxd**2 + vyd**2 + vzd**2)
        
        # compute relative position vectors / baselines
        RelPos[:,:,i,j] = np.concatenate((xd[:,0].reshape(L,1) - xd[:,1].reshape(L,1) ,
                                        yd[:,0].reshape(L,1) - yd[:,1].reshape(L,1),
                                        zd[:,0].reshape(L,1) - zd[:,1].reshape(L,1)),
                                       axis=1)
         
        if i > 0:
            Dmoon[:,:,i,j] = moon[:,:,0,j] - moon[:,:,i,j]
            Dx[:,:,i,j] = x[:,:,0,j] - x[:,:,i,j]
            Dy[:,:,i,j] = y[:,:,0,j] - y[:,:,i,j]
            Dz[:,:,i,j] = z[:,:,0,j] - z[:,:,i,j]
            Dvx[:,:,i,j] = vx[:,:,0,j] - vx[:,:,i,j]
            Dvy[:,:,i,j] = vy[:,:,0,j] - vy[:,:,i,j]
            Dvz[:,:,i,j] = vz[:,:,0,j] - vz[:,:,i,j]
            
            Dr[:,:,i,j] = np.sqrt(Dx[:,:,i,j]**2 + Dy[:,:,i,j]**2 + Dz[:,:,i,j]**2)
            Dv[:,:,i,j] = np.sqrt(Dvx[:,:,i,j]**2 + Dvy[:,:,i,j]**2 + Dvz[:,:,i,j]**2)
            
            
            DRelPos[:,:,i,j] = RelPos[:,:,0,j] - RelPos[:,:,i,j]
            MDRelPos[:,0,i,j] = np.sqrt(DRelPos[:,0,i,j]**2 + DRelPos[:,1,i,j]**2 + DRelPos[:,2,i,j]**2)
            
            MDmoon[:,0,i,j] = np.sqrt(Dmoon[:,0,i,j]**2 + Dmoon[:,1,i,j]**2 + Dmoon[:,2,i,j]**2)
        
print("Succesfully loaded all data")
del xd, yd, zd, vxd, vyd, vzd

t -= t[0]

#%%  Plot preparation
Dr[Dr == 0] = np.NaN

PF.ClearPlots()

legend = ['Default model',
          'Mars point mass removed',
          'Venus point mass removed',
          'Solar PM',
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
# d = LL.folder_dict_funcs["PA"]

# For performance reasons subdivide the data before displaying the wireframe
subdiv = 300 
# fig,ax = PF.MakeFig((14,8), '3D')
xp = t * np.ones((6,L))
yp = [1e3, 10e3, 30e3, 50e3, 70e3, 100e3]* np.ones((L,6))
zp = MDRelPos[:,0,3,:] # change in x of baselines
# PF.Simple3DWireframe(xp.T[::subdiv], yp[::subdiv], zp[::subdiv], xmod = 24*3600., ymod = 1000., zmod = 1000.,
#                  fig = fig, ax = ax)
# PF.DressFig(fig, ax, title="Total baseline magnitude change from removal of solar point mass", xlabel = "Time [days",
#             ylabel = "Initial x-displacement [km]", zlabel="Baseline magnitude change [km]")

# PF.Animate3DPlotToGIF(fig, ax, "testgif", d)

#%%
# fig,ax = PF.MakeFig((14,8), '3D')
# zp = DRelPos[:,0,3,:] # change in x of baselines
# PF.Simple3DWireframe(xp.T[::subdiv], yp[::subdiv], zp[::subdiv], xmod = 24*3600., ymod = 1000., zmod = 1000.,
#                  fig = fig, ax = ax)
# PF.DressFig(fig, ax, title="Total baseline x direction change from removal of solar point mass", xlabel = "Time [days",
#             ylabel = "Initial x-displacement [km]", zlabel="Change in baseline x component [km]")


# fig,ax = PF.MakeFig((14,8), '3D')
# zp = DRelPos[:,1,3,:] # change in y of baselines
# PF.Simple3DWireframe(xp.T[::subdiv], yp[::subdiv], zp[::subdiv], xmod = 24*3600., ymod = 1000., zmod = 1000.,
#                  fig = fig, ax = ax)
# PF.DressFig(fig, ax, title="Total baseline y direction change from removal of solar point mass", xlabel = "Time [days",
#             ylabel = "Initial x-displacement [km]", zlabel="Change in baseline y component [km]")

# fig,ax = PF.MakeFig((14,8), '3D')
# zp = DRelPos[:,0,3,:] # change in z of baselines
# PF.Simple3DWireframe(xp.T[::subdiv], yp[::subdiv], zp[::subdiv], xmod = 24*3600., ymod = 1000., zmod = 1000.,
#                  fig = fig, ax = ax)
# PF.DressFig(fig, ax, title="Total baseline z direction change from removal of solar point mass", xlabel = "Time [days",
#             ylabel = "Initial x-displacement [km]", zlabel="Change in baseline z component [km]")

#%% Linear regression fit
dataset = 3

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(MDRelPos[-1,0,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()
print(fittedmodel.params)
print(fittedmodel.summary() )

fig,ax = PF.MakeFig()

Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax,xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax)
PF.DressFig(fig, ax, title="Baseline magnitude change from " + legend[dataset] + " after 1 year",
            xlabel='Initial baseline [km]', ylabel='Magnitude change in baseline after 1 year [km]')
PF.saveFig(d, "R linearity")

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(DRelPos[-1,0,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()
print(fittedmodel.params)
print(fittedmodel.summary() )

fig,ax = PF.MakeFig()

Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax,xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax)
PF.DressFig(fig, ax, title="Baseline X change from " + legend[dataset] + " after 1 year",
            xlabel='Initial baseline [km]', ylabel='X coordinate change [km]')
PF.saveFig(d, "X linearity")

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(DRelPos[-1,1,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()
print(fittedmodel.params)
print(fittedmodel.summary() )

fig,ax = PF.MakeFig()

Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax,xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax)
PF.DressFig(fig, ax, title="Baseline Y change from " + legend[dataset] + " after 1 year",
            xlabel='Initial baseline [km]', ylabel='Y coordinate change [km]')
PF.saveFig(d, "Y linearity")

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(DRelPos[-1,2,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()
print(fittedmodel.params)
print(fittedmodel.summary() )

fig,ax = PF.MakeFig()

Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax,xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax)
PF.DressFig(fig, ax, title="Baseline Z change from " + legend[dataset] + " after 1 year",
            xlabel='Initial baseline [km]', ylabel='Z coordinate change [km]')
PF.saveFig(d, "Z linearity")

#%% Compound lin. regression fit in 1 fig

dataset = 3

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(MDRelPos[-1,0,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()
print(fittedmodel.params)
print(fittedmodel.summary() )

fig,ax = pp.subplots(1,3, figsize = (12,4))
# pp.suptitle("Change in baseline coordinates from Solar PM after 1 year in orbit")

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(DRelPos[-1,0,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()


Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax[0],xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax[0])
PF.DressFig(fig, ax[0], title="",
            xlabel='Initial baseline [km]', ylabel='X coordinate change [km]')
# # PF.saveFig(d, "X linearity")
# R squared:{np.round(fittedmodel.rsquared,8)}
X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(DRelPos[-1,1,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()

Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax[1],xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax[1])
PF.DressFig(fig, ax[1], title='',
            xlabel='Initial baseline [km]', ylabel='Y coordinate change [km]')
# PF.saveFig(d, "Y linearity")

X = yp[0,:]
X = sm.add_constant(X)
Y = []
for i in range(0,len(addendumList)):
    Y.append(DRelPos[-1,2,dataset,i])
    
lin_reg = sm.OLS(Y,X)
fittedmodel = lin_reg.fit()

Y = np.array(Y)
PF.Simple2DScatter(X[:,1], Y, fig = fig, ax = ax[2],xmod = 1000, ymod = 1000)
modelplot = fittedmodel.params[0] + fittedmodel.params[1]*X[:,1]
PF.Simple2DPlot(X[:,1],modelplot,xmod = 1000, ymod = 1000, fig=fig,ax = ax[2])
PF.DressFig(fig, ax[2], title="",
            xlabel='Initial baseline [km]', ylabel='Z coordinate change [km]')

fig.tight_layout()

pp.savefig(d + "Linearity_compound")
# # PF.saveFig(d, "Z linearity")
