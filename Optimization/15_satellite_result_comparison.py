# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 15:59:02 2020

@author: mrfan
"""

import numpy as np
import sys 

from os import chdir
from os import getcwd

from matplotlib import pyplot as pp
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

sys.path.append('../')
chdir('../')
import LazyLib as LL
import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF



#%%

figure_subfolder = f"/15sats365/"
figure_folder = LL.folder_dict["Results"] + figure_subfolder

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization/Coarse_Multi_365/'
#%% Import data
PF.CloseAll()
gen = 76

made = False
data = []   
for i in range(0,7):
    
    fitness = f'fitness_singlemethod_long_sd42_sats15_nisl8_npop32_int1_g{gen}_i{i}.dat'
    population = f'population_singlemethod_long_sd42_sats15_nisl8_npop32_int1_g{gen}_i{i}.dat'
    
    fitness = DI.ImportSingleColumnFile(datafolder + fitness)
    popdata = DI.ImportDataFile(datafolder + population)
    
    for j in range(0,len(fitness)):
        if fitness[j] == 0: 
            if made:
                data = np.append(data, popdata[j,:].reshape(1,len(popdata[0,:])) ,axis = 0 )  
            else: 
                data = np.array(popdata[j,:]).reshape(1,len(popdata[0,:]))
                made = True

t,x,y,z,vx,vy,vz,ddd = DI.ImportPropagationHistory(datafolder + 
                                             "propagationHistory_singlemethod_long_sd42_sats15_nisl8_npop32_int1_best.dat")

moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")

sat = np.concatenate((x,y,z), axis=1)
x_B,y_B,z_B, L4B = DI.ConvertToBarycentricFrame(moon, sat,exportL4 = True)

core = DI.ImportDataFile(datafolder + 
                         "corePosition_singlemethod_long_sd42_sats15_nisl8_npop32_int1_best.dat")

core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core[:3].flatten())
core_x_Bf , core_y_Bf, core_z_Bf, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)

#%% Extract and convert data to barycentric relative frame    
xdata = data[:,6::3]
ydata = data[:,7::3]
zdata = data[:,8::3]

L = len(data[:,0])
moon = moon[0,:]*np.ones((L*15,6))
L4 = DI.computeL4Location(moon)


core_x_J2000 = (data[:,0] + L4[0,0]).reshape(L,1)*np.ones((1,15))
core_y_J2000 = (data[:,1] + L4[0,1]).reshape(L,1)*np.ones((1,15))
core_z_J2000 = (data[:,2] + L4[0,2]).reshape(L,1)*np.ones((1,15))
xdata_J2000 = xdata + core_x_J2000
ydata_J2000 = ydata + core_y_J2000
zdata_J2000 = zdata + core_z_J2000

L = len(core_x_J2000.flatten())
cores_J2000 = np.zeros((L,3))
cores_J2000[:,0] = core_x_J2000.flatten()
cores_J2000[:,1] = core_y_J2000.flatten()
cores_J2000[:,2] = core_z_J2000.flatten()


sats_J2000 = np.zeros((L,3))
sats_J2000[:,0] = xdata_J2000.flatten()
sats_J2000[:,1] = ydata_J2000.flatten()
sats_J2000[:,2] = zdata_J2000.flatten()

core_x_B, core_y_B, core_z_B,d = DI.ConvertToBarycentricFrame(moon, cores_J2000)
xdata_B, ydata_B, zdata_B,d = DI.ConvertToBarycentricFrame(moon, sats_J2000)
L4_x_B, L4_y_B, L4_z_B,d = DI.ConvertToBarycentricFrame(moon, L4)

rel_x_B = xdata_B - core_x_B
rel_y_B = ydata_B - core_y_B
rel_z_B = zdata_B - core_z_B

rel_cx_B = core_x_B[::15] - L4_x_B[::15]
rel_cy_B = core_y_B[::15] - L4_y_B[::15]
rel_cz_B = core_z_B[::15] - L4_z_B[::15]
#%% Convert core velocities to barycentric 

ang = np.cross(moon[0,:3],moon[0,3:])/(np.linalg.norm(moon[0,:3])**2)
VL4 = np.cross(ang,L4[0])

rel_vx_J2000 = data[:,3]
rel_vy_J2000 = data[:,4]
rel_vz_J2000 = data[:,5]

H  = np.cross(moon[0,:3],moon[0,3:])
b_x = moon[0,:3]/np.linalg.norm(moon[0,:3])
b_z = H/np.linalg.norm(H)
b_y = np.cross(b_z,b_x)

rel_vx_B, rel_vy_B, rel_vz_B = (np.empty_like(rel_vx_J2000) for i in range(0,3))
b_vel = np.zeros((len(data[:,3]),3))
for j in range(0,len(data[:,0])):
    
    rel_v = np.array([rel_vx_J2000[j], rel_vy_J2000[j], rel_vz_J2000[j]])
    
    rel_vx_B[j] = GF.vectorProj(rel_v,b_x)
    rel_vy_B[j] = GF.vectorProj(rel_v,b_y)
    rel_vz_B[j] = GF.vectorProj(rel_v,b_z)
    
    b_vel[j,:] = np.array([rel_vx_B[i], rel_vy_B[i], rel_vz_B[i]])


#%% Analyse data clusters
                
fig,ax = PF.MakeFig(dimensions = '3D')
ax.scatter(data[:,0]/1e3,data[:,1]/1e3,data[:,2]/1e3, color='r')
ax.quiver(data[:,0]/1e3,data[:,1]/1e3,data[:,2]/1e3, data[:,3],data[:,4],data[:,5],length= 0.2)
ax.scatter(0,0,0,s = 100, marker = 's')

#%%

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.scatter(data[:,0]/1e3,data[:,1]/1e3,marker='s', color = 'r')
ax1.scatter(0,0,s = 100, marker = 'd', color = 'g')
ax1.quiver(data[:,0]/1e3,data[:,1]/1e3, data[:,3],data[:,4])

# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('x [km]')
ax1.set_ylabel('y [km]')

ax3.scatter(data[:,1]/1e3,data[:,2]/1e3,marker='s', color = 'r')
ax3.scatter(0,0,s = 100, marker = 'd', color = 'g')
ax3.quiver(data[:,1]/1e3,data[:,2]/1e3, data[:,4],data[:,5])
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('y [km]')
ax3.set_ylabel('z [km]')

ax2.scatter(data[:,0]/1e3,data[:,2]/1e3,marker='s', color = 'r')
ax2.scatter(0,0,s = 100, marker = 'd', color = 'g')
ax2.quiver(data[:,0]/1e3,data[:,2]/1e3, data[:,3],data[:,5])
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('x [km]')
ax2.set_ylabel('z [km]')

fig.legend(['Core positions', 'L4', 'Core velocity' ],bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_solution_distributions")
#%% Same figure but for barycentric frame


fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.scatter(rel_x_B/1e3,rel_y_B/1e3, color='g', alpha = 0.2,zorder=9)
ax1.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax1.quiver(np.zeros(rel_vy_B.shape),np.zeros(rel_vy_B.shape),
            rel_vx_B,rel_vy_B,scale = 60,alpha=0.3,color='k')

# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter(rel_y_B/1e3,rel_z_B/1e3, color='g', alpha = 0.2,zorder=9)
ax3.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)

ax3.quiver(np.zeros(rel_vy_B.shape),np.zeros(rel_vz_B.shape),
            rel_vy_B,rel_vz_B,scale = 60,alpha=0.3,color='k')
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter(rel_x_B/1e3,rel_z_B//1e3, color='g', alpha = 0.2,zorder=9)
ax2.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax2.quiver(np.zeros(rel_vx_B.shape),np.zeros(rel_vz_B.shape),
            rel_vx_B,rel_vz_B,scale = 60,alpha=0.3,color='k')
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')

fig.legend(['Satellite positions', 'Core', 'Core velocity' ],bbox_to_anchor=(0.52, 0.97))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_satellite_dist_Barycentric")

#%% Barycentric core distribution

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.scatter(rel_cx_B/1e3,rel_cy_B/1e3, color='r', alpha = 0.8,zorder=9,marker='s')
ax1.scatter(0,0,s = 100, marker = 'd', color = 'b',zorder=10)
ax1.quiver(rel_cx_B/1e3,rel_cy_B/1e3,
            rel_vx_B,rel_vy_B,scale = 60,alpha=0.3,color='k')

# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter(rel_cy_B/1e3,rel_cz_B/1e3, color='r', alpha = 0.8,zorder=9,marker='s')
ax3.scatter(0,0,s = 100, marker = 'd', color = 'b',zorder=10)

ax3.quiver(rel_cy_B/1e3,rel_cz_B/1e3,
            rel_vy_B,rel_vz_B,scale = 60,alpha=0.3,color='k')
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter(rel_cx_B/1e3,rel_cz_B/1e3, color='r', alpha = 0.8,zorder=9,marker='s')
ax2.scatter(0,0,s = 100, marker = 'd', color = 'b',zorder=10)
ax2.quiver(rel_cx_B/1e3,rel_cz_B/1e3,
            rel_vx_B,rel_vz_B,scale = 60,alpha=0.3,color='k')
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')

fig.legend(['Core positions', 'L4', 'Core velocity' ],bbox_to_anchor=(0.51, 0.97))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_core_dist_Barycentric")

#%% J2000 swarm distribution


fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.scatter(xdata.flatten()/1e3,ydata.flatten()/1e3, color='g', alpha = 0.2,zorder=9)
ax1.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax1.quiver(np.zeros(data[:,4].shape),np.zeros(data[:,4].shape),
           data[:,3],data[:,4],scale = 60,alpha=0.3,color='k')

# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter(ydata.flatten()/1e3,zdata.flatten()/1e3, color='g', alpha = 0.2,zorder=9)
ax3.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)

ax3.quiver(np.zeros(data[:,4].shape),np.zeros(data[:,4].shape),
           data[:,4],data[:,5],scale = 60,alpha=0.3,color='k')
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter(xdata.flatten()/1e3,zdata.flatten()/1e3, color='g', alpha = 0.2,zorder=9)
ax2.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax2.quiver(np.zeros(data[:,3].shape),np.zeros(data[:,4].shape),
           data[:,3],data[:,5],scale = 60,alpha=0.3,color='k')
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')

fig.legend(['Satellite positions', 'Core', 'Core velocity' ],bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_satellite_dist_J2000")

#%% 3d relative satellite plot

fig,ax = PF.MakeFig(dimensions='3D')

ax.scatter(0,0,0,marker='s',s = 120)
ax.scatter(rel_x_B/1e3, rel_y_B/1e3, rel_z_B/1e3, alpha = 0.5)
ax.scatter(rel_x_B[:15]/1e3, rel_y_B[:15]/1e3, rel_z_B[:15]/1e3,color='k', alpha = 1)


#%% 3d orbit plot

fig,ax = PF.MakeFig(figsize = (12,4) , dimensions='3D')

ax.scatter(0,0,0,marker='d',s = 120,color='b')

ax.plot((x_B[:,0]-L4B[:,0])/1e3,(y_B[:,0]-L4B[:,1])/1e3,(z_B[:,0]-L4B[:,2])/1e3)
ax.scatter((x_B[-1,0]-L4B[-1,0])/1e3,(y_B[-1,0]-L4B[-1,1])/1e3,
           (z_B[-1,0]-L4B[-1,2])/1e3, marker='s',s = 120,color='r')
ax.set_xlabel('$x_B$ [km]')
ax.set_ylabel('$y_B$ [km]')
ax.set_zlabel('$z_B$ [km]')
pp.legend(["L4", "Swarm motion", "Final position"])

pp.savefig(figure_folder + "3d-swarm-motion")


#%% Initial motion rel. to core

cutoff = 1*30*24
core = core[:3].flatten()

rel_x_B = (x_B - core_x_Bf)/1e3
rel_y_B = (y_B - core_y_Bf)/1e3
rel_z_B = (z_B - core_z_Bf)/1e3


fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))




ax1.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10,label= 'Core position')
ax1.scatter(rel_x_B[0,:],rel_y_B[0,:], color='k', alpha = 1,zorder=9, label = 'Init. satellite position')

ax1.plot(rel_x_B[:cutoff,0], rel_y_B[:cutoff,0], color = 'b', alpha = 0.6,label='Satellite motion 30d')
ax3.plot(rel_y_B[:cutoff,0], rel_z_B[:cutoff,0], color = 'b', alpha = 0.6)
ax2.plot(rel_x_B[:cutoff,0], rel_z_B[:cutoff,0], color = 'b', alpha = 0.6)
for i in range(1,15):
    ax1.plot(rel_x_B[:cutoff,i], rel_y_B[:cutoff,i], color = 'b', alpha = 0.6)
    ax3.plot(rel_y_B[:cutoff,i], rel_z_B[:cutoff,i], color = 'b', alpha = 0.6)
    ax2.plot(rel_x_B[:cutoff,i], rel_z_B[:cutoff,i], color = 'b', alpha = 0.6)
# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)

ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax3.scatter(rel_y_B[0,:],rel_z_B[0,:], color='k', alpha = 1,zorder=9)
# ax3.quiver(np.zeros(data[:,4].shape),np.zeros(data[:,4].shape),
#            data[:,4],data[:,5],scale = 60,alpha=0.3,color='k')
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax2.scatter(rel_x_B[0,:],rel_z_B[0,:], color='k', alpha = 1,zorder=9)
# ax2.quiver(np.zeros(data[:,3].shape),np.zeros(data[:,4].shape),
#            data[:,3],data[:,5],scale = 60,alpha=0.3,color='k')
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')


fig.legend(bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_satellite_rel_motion_B_1m")

#%%


fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))
ax1.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10,label= 'Core position')
ax1.scatter(rel_x_B[0,:],rel_y_B[0,:], color='k', alpha = 1,zorder=11, label = 'Init. satellite position')

ax1.plot(rel_x_B[:2160,0], rel_y_B[:2160,0], color = 'b', alpha = 0.6,label='Satellite motion 90d')
ax3.plot(rel_y_B[:2160,0], rel_z_B[:2160,0], color = 'b', alpha = 0.6)
ax2.plot(rel_x_B[:2160,0], rel_z_B[:2160,0], color = 'b', alpha = 0.6)
ax1.plot(rel_x_B[:cutoff,0], rel_y_B[:cutoff,0], color = 'y', alpha = 1,zorder=9,label='Satellite motion 30d')
ax3.plot(rel_y_B[:cutoff,0], rel_z_B[:cutoff,0], color = 'y', alpha = 1,zorder=9)
ax2.plot(rel_x_B[:cutoff,0], rel_z_B[:cutoff,0], color = 'y', alpha = 1,zorder=9)
for i in range(1,15):
    ax1.plot(rel_x_B[:2160,i], rel_y_B[:2160,i], color = 'b', alpha = 0.6)
    ax3.plot(rel_y_B[:2160,i], rel_z_B[:2160,i], color = 'b', alpha = 0.6)
    ax2.plot(rel_x_B[:2160,i], rel_z_B[:2160,i], color = 'b', alpha = 0.6)
    ax1.plot(rel_x_B[:cutoff,i], rel_y_B[:cutoff,i], color = 'y', alpha = 1,zorder=9)
    ax3.plot(rel_y_B[:cutoff,i], rel_z_B[:cutoff,i], color = 'y', alpha = 1,zorder=9)
    ax2.plot(rel_x_B[:cutoff,i], rel_z_B[:cutoff,i], color = 'y', alpha = 1,zorder=9)
# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)

ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax3.scatter(rel_y_B[0,:],rel_z_B[0,:], color='k', alpha = 1,zorder=11)
# ax3.quiver(np.zeros(data[:,4].shape),np.zeros(data[:,4].shape),
#            data[:,4],data[:,5],scale = 60,alpha=0.3,color='k')
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax2.scatter(rel_x_B[0,:],rel_z_B[0,:], color='k', alpha = 1,zorder=11)
# ax2.quiver(np.zeros(data[:,3].shape),np.zeros(data[:,4].shape),
#            data[:,3],data[:,5],scale = 60,alpha=0.3,color='k')
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')


fig.legend(bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_satellite_rel_motion_B_3m")


#%%

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))
ax1.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10,label= 'Core position')
ax1.scatter(rel_x_B[0,:],rel_y_B[0,:], color='k', alpha = 1,zorder=11, label = 'Init. satellite position')

ax1.plot(rel_x_B[:,0], rel_y_B[:,0], color = 'b', alpha = 0.6,label='Satellite motion 365d')
ax3.plot(rel_y_B[:,0], rel_z_B[:,0], color = 'b', alpha = 0.6)
ax2.plot(rel_x_B[:,0], rel_z_B[:,0], color = 'b', alpha = 0.6)
ax1.plot(rel_x_B[:cutoff,0], rel_y_B[:cutoff,0], color = 'y', alpha = 1,zorder=9,label='Satellite motion 30d')
ax3.plot(rel_y_B[:cutoff,0], rel_z_B[:cutoff,0], color = 'y', alpha = 1,zorder=9)
ax2.plot(rel_x_B[:cutoff,0], rel_z_B[:cutoff,0], color = 'y', alpha = 1,zorder=9)
for i in range(1,15):
    ax1.plot(rel_x_B[:,i], rel_y_B[:,i], color = 'b', alpha = 0.6)
    ax3.plot(rel_y_B[:,i], rel_z_B[:,i], color = 'b', alpha = 0.6)
    ax2.plot(rel_x_B[:,i], rel_z_B[:,i], color = 'b', alpha = 0.6)
    ax1.plot(rel_x_B[:cutoff,i], rel_y_B[:cutoff,i], color = 'y', alpha = 1,zorder=9)
    ax3.plot(rel_y_B[:cutoff,i], rel_z_B[:cutoff,i], color = 'y', alpha = 1,zorder=9)
    ax2.plot(rel_x_B[:cutoff,i], rel_z_B[:cutoff,i], color = 'y', alpha = 1,zorder=9)
# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)

ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax3.scatter(rel_y_B[0,:],rel_z_B[0,:], color='k', alpha = 1,zorder=11)
# ax3.quiver(np.zeros(data[:,4].shape),np.zeros(data[:,4].shape),
#            data[:,4],data[:,5],scale = 60,alpha=0.3,color='k')
# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter(0,0,s = 100, marker = 's', color = 'r',zorder=10)
ax2.scatter(rel_x_B[0,:],rel_z_B[0,:], color='k', alpha = 1,zorder=11)
# ax2.quiver(np.zeros(data[:,3].shape),np.zeros(data[:,4].shape),
#            data[:,3],data[:,5],scale = 60,alpha=0.3,color='k')
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')


fig.legend(bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_satellite_rel_motion_B_1y")

#%%
td = (t-t[0])/(24*3600)
fig, (ax1,ax2,ax3) = pp.subplots(3,1, figsize = (15,7))

for i in range(0,15):
    ax1.plot(td,rel_x_B[:,i],color='b', alpha = 0.5)
    ax2.plot(td,rel_y_B[:,i],color='b', alpha = 0.5)
    ax3.plot(td,rel_z_B[:,i],color='b', alpha = 0.5)

ax1.set_ylabel('Relative $x_B$[km]')
ax2.set_ylabel('Relative $y_B$[km]')

ax3.set_xlabel('time in orbit [days]')
ax3.set_ylabel('Relative $z_B$[km]')

ax1.set_xlim([0,365])
ax2.set_xlim([0,365])
ax3.set_xlim([0,365])

pp.tight_layout()
pp.savefig(figure_folder + 'harmonic_decay')