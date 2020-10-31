# -*- coding: utf-8 -*-
"""
Compare trajectories between body-pointing and butt-pointing satellites
Created on Wed Mar 25 13:26:46 2020

@author: mrfan
"""
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

sys.path.append('../')
chdir('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF
from matplotlib import pyplot as pp
from matplotlib.ticker import FormatStrFormatter
chdir('PerturbationAnalysis/')

RE = LL.Constants["Lunar_mean_dist"];
mu = LL.Constants["mu*"]

def combinedOrbitPlot(core_x,core_y,core_z,filename):
    xlim, ylim, zlim = [-.4, 1.2], [-.1,1.3], [-2, 2]
    
    fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

    
    ax2.plot(core_x/RE,core_y/RE)
    ax2.scatter(-mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
    ax2.scatter(1-mu,0,marker = 'o',c = 'g', zorder=10)
    ax2.scatter(0.5,np.sqrt(3)/2,marker = 'd', c = 'k', zorder=10)
    ax2.set_xlabel('x [R]')
    ax2.set_ylabel('y [R]')
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)
    
    
    ax1.plot(core_y/RE,core_z*100/RE)
    ax1.scatter(0,0,marker = 'o', s = 70, c = 'r', zorder=10)
    ax1.scatter(0,0,marker = 'o',c = 'g', zorder=10)
    ax1.scatter(np.sqrt(3)/2,0,marker = 'd', c = 'k', zorder=10)
    ax1.set_xlabel('y [R]')
    ax1.set_ylabel('z [R * 1E-2]')
    ax1.set_xlim(ylim)
    ax1.set_ylim(zlim)
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    
    ax3.plot(core_x/RE,core_z*100/RE)
    ax3.scatter(-mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
    ax3.scatter(1-mu,0,marker = 'o',c = 'g', zorder=10)
    ax3.scatter(0.5,0,marker = 'd', c = 'k', zorder=10)
    ax3.set_xlabel('x [R]')
    ax3.set_ylabel('z [R * 1E-2]')
    ax3.set_xlim(xlim)
    ax3.set_ylim(zlim)
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    fig.legend(['Swarm \n motion', 'Earth', 'Moon', 'L4' ],loc=(0.905,0.5))
    
    pp.savefig(d + filename)
    return


#%% SETTINGS
# Draw 3d plots of motion relative to L4 in barycentric frame, needs large time scales to show effects
barycentricL4MotionComparison = True;

addendum = "_180d"
# Number of cases to be taken into account
numcases = 2
# Moon-model affecting border
moonborder = 13


#%% Data import
chdir("../")
cwd = getcwd()
append = cwd + "/Data/"

print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")

# x,y,z,vx,vy,vz, R, V = (np.ndarray([]) for i in range(0,8))

datafolder = cwd + '/Data/PointingAnalysis/'

moon = np.ndarray([])

madevars = False
for i in range( 0 , numcases):
    filename = datafolder + "propagationHistory_" + str(i) + ".dat"
    
    t, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,0, False)
    md = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")
    if madevars == False:
        moon = np.zeros((md.shape[0],md.shape[1],numcases))
        Dmoon = np.zeros((md.shape[0],md.shape[1],numcases))
        x,y,z,vx,vy,vz, R,V, E  = (np.zeros((xd.shape[0],xd.shape[1],numcases)) 
                                   for i in range(0,9))
        Dx,Dy,Dz,Dvx,Dvy,Dvz,Dr,Dv,DE  = (np.zeros((xd.shape[0],xd.shape[1],numcases)) 
                                       for i in range(0,9))
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
        
print("Succesfully loaded all data")
del xd, yd, zd, vxd, vyd, vzd

#%%  Plot preparation

PF.ClearPlots()

legend = ['Drift in x-direction',
          'Drift in y-direction',
          'Drift in z-direction',
          'Total positional drift'
    ]

styles = [None, None, '-d', None,'-o', '--', '--', '--', '--', '-.','-.','-.', None, None, None]

#%% Do the plotting
d = LL.folder_dict["Results"]
#d = getcwd() + '/Figures/Perturbation Analysis/'
LL.CheckDir(d)

fig, ax = PF.Matrix2DPlot(t,Dx[:,1,1].reshape(Dx.shape[0], 1, 1),xmod=24*3600, figsize=(12,4))
PF.Matrix2DPlot(t,Dy[:,1,1].reshape(Dx.shape[0], 1, 1),xmod=24*3600, fig = fig, ax = ax)
PF.Matrix2DPlot(t,Dz[:,1,1].reshape(Dx.shape[0], 1, 1),xmod=24*3600, fig = fig, ax = ax)
PF.Matrix2DPlot(t,Dr[:,1,1].reshape(Dx.shape[0], 1, 1), title='Changes in position due to body pointing',
                xlabel = 't[days]',ylabel='change [km]', xmod=24*3600,
                fig=fig, ax = ax,
                savefig = True, figfolder = d, name='Changes_Satellite_Modelling',
                legendlist = legend, stylelist = styles, markN = 50,xlim = [0,365],
                logScale = False)

fig, ax = PF.Matrix2DPlot(t,np.abs(Dx[:,1,1].reshape(Dx.shape[0], 1, 1)),xmod=24*3600, logScale = True, figsize=(12,4))
PF.Matrix2DPlot(t,np.abs(Dy[:,1,1].reshape(Dx.shape[0], 1, 1)),xmod=24*3600, fig = fig, ax = ax)
PF.Matrix2DPlot(t,np.abs(Dz[:,1,1].reshape(Dx.shape[0], 1, 1)),xmod=24*3600, fig = fig, ax = ax)
PF.Matrix2DPlot(t,np.abs(Dr[:,1,1].reshape(Dx.shape[0], 1, 1)), title='Changes in absolute position due to body pointing',
                xlabel = 't[days]',ylabel='change [km]', xmod=24*3600,
                fig=fig, ax = ax,
                savefig = True, figfolder = d, name='Changes_Satellite_Modelling_Absolute',
                legendlist = legend, stylelist = styles, markN = 50,xlim = [0,365],
                logScale = True)

PF.Matrix2DPlot(t,DE[:,1:,1:], title='(Absolute)Changes in system energy',
                xlabel = 't[days]',ylabel='$\Delta$E [kJ]',xmod=24*3600,
                savefig = True, figfolder = d, name='Changes_System_Energy',xlim = [0,365],
                legendlist = ["Change in satellite potential energy from body pointing"], stylelist = styles,markN = 50,
                logScale = False, figsize=(12,4))
#%%
if barycentricL4MotionComparison:
    for i in range(0,numcases):
        sat = np.concatenate([x[:,:,i],y[:,:,i],z[:,:,i]],axis=1)
        co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon[:,:,i],sat, exportL4 =True)
        combinedOrbitPlot(co_x[:,1],co_y[:,1],co_z[:,1],f"combined_orbit_S2_{i}")
        N = co_x.shape[0]
        fig2, ax2 = PF.Simple3DPlot(co_x - L4b[:,0].reshape(N,1)*np.ones(co_x.shape), 
                                    co_y - L4b[:,1].reshape(N,1)*np.ones(co_y.shape),
                                    co_z - L4b[:,2].reshape(N,1)*np.ones(co_z.shape))
        PF.Simple3DScatter(0,0,0, fig= fig2, ax = ax2)
            
#%%
        
PF.Simple3DPlot(Dx[:,1,1], Dy[:,1,1], Dz[:,1,1], xlabel='$\Delta$ x [km]', ylabel = '$\Delta$ y [km]',
                zlabel = '$\Delta$ z [km]')

PF.saveFig(d, "3d-position-change")
#%%


fig, ax = PF.Matrix2DPlot(t[:480],np.abs(Dx[:480,1,1].reshape(480, 1, 1)),xmod=3600, logScale = True, figsize=(12,4))
PF.Matrix2DPlot(t[:480],np.abs(Dy[:480,1,1].reshape(480, 1, 1)),xmod=3600, fig = fig, ax = ax)
PF.Matrix2DPlot(t[:480],np.abs(Dz[:480,1,1].reshape(480, 1, 1)),xmod=3600, fig = fig, ax = ax)

pp.xticks(np.arange(t[0]/3600, t[480]/3600, 12.0))
PF.Matrix2DPlot(t[:480],np.abs(Dr[:480,1,1].reshape(480, 1, 1)), title='Changes in absolute position due to body pointing',
                xlabel = 't [hours]',ylabel='change [km]', xmod=3600,
                fig=fig, ax = ax,
                savefig = True, figfolder = d, name='Changes_Satellite_Modelling_Absolute_10d',
                legendlist = legend, stylelist = styles, markN = 50,xlim = [0,240],
                logScale = True)