# -*- cofding: utf-8 -*-
"""
Detailed single-solution analysis script, for the real-deal long term propagated orbits
Created on Wed Jul 22 15:07:57 2020

@author: Jurriaan van 't Hoff'
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


def combinedOrbitPlot(core_x,core_y,core_z,filename, cof= 0):
    xlim, ylim, zlim = [-.4, 1.2], [-.1,1.3], [-2, 2]
    
    fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))
    
    ax1.plot(core_x/RE,core_y/RE)
    if cof != 0: ax1.plot(core_x[:cof]/RE,core_y[:cof]/RE,color= 'y')
    ax1.scatter(-mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
    ax1.scatter(1-mu,0,marker = 'o',c = 'g', zorder=10)
    ax1.scatter(0.5,np.sqrt(3)/2,marker = 'd', c = 'k', zorder=10)
    ax1.set_xlabel('$x_B$ [R]')
    ax1.set_ylabel('$y_B$ [R]')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    
    
    ax3.plot(core_y/RE,core_z*100/RE)
    if cof != 0: ax3.plot(core_y[:cof]/RE,core_z[:cof]*100/RE,color= 'y')
    ax3.scatter(0,0,marker = 'o', s = 70, c = 'r', zorder=10)
    ax3.scatter(0,0,marker = 'o',c = 'g', zorder=10)
    ax3.scatter(np.sqrt(3)/2,0,marker = 'd', c = 'k', zorder=10)
    ax3.set_xlabel('$y_B$ [R]')
    ax3.set_ylabel('$z_B$ [R * 1E-2]')
    ax3.set_xlim(ylim)
    ax3.set_ylim(zlim)
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    
    ax2.plot(core_x/RE,core_z*100/RE)
    if cof != 0: ax2.plot(core_x[:cof]/RE,core_z[:cof]*100/RE,color= 'y')
    ax2.scatter(-mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
    ax2.scatter(1-mu,0,marker = 'o',c = 'g', zorder=10)
    ax2.scatter(0.5,0,marker = 'd', c = 'k', zorder=10)
    ax2.set_xlabel('$x_B$ [R]')
    ax2.set_ylabel('$z_B$ [R * 1E-2]')
    ax2.set_xlim(xlim)
    ax2.set_ylim(zlim)
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    if cof !=0:
        fig.legend(['Swarm \n motion 5yr', 'Swarm \n motion 1yr', 'Earth', 'Moon', 'L4' ],loc=(0.71,0.64))
    else:
        fig.legend(['Swarm \n motion', 'Earth', 'Moon', 'L4' ],loc=(0.71,0.7))
    pp.tight_layout()
    pp.savefig(figure_folder + filename)
    return

def vectorProj(v1,v2):
    """Projection of vector v1 onto v2"""
    
    v3 = np.dot(v1,v2)/np.linalg.norm(v2)
    return v3

Nsats = 35

figure_subfolder = f"/{Nsats}SatsChampion_2nd_25A2/"
data_subfolder = "/champions_propagated/35sats2ndorder25A2/" #35sats2ndorder25A2/

set_name = f"{Nsats}sat_champ"


figure_folder = LL.folder_dict["Results"] + figure_subfolder

LL.CheckDir(figure_folder)

cwd = getcwd()
datafolder = cwd + '/Data/Optimization' + data_subfolder

#%% 
PF.CloseAll()
display_cutoff = 365

cof = int(display_cutoff*24/4)

file_history = f"propagationHistory_{set_name}.dat"
file_core = f"corePosition_{set_name}.dat"
# compute initial satellite positions relative to core
t, x, y, z, vx, vy, vz,ddd  = DI.ImportPropagationHistory(datafolder + file_history,0, False)
sat = np.concatenate((x,y,z), axis=1)


moon = DI.ImportMoonHistory(datafolder + f"propagationHistory_{set_name}_moon.dat")
L4 = DI.computeL4Location(moon)

core = DI.ImportSingleColumnFile(datafolder + file_core)
core_pos = GF.ExtrapolateSwarmMotion(x, y, z, core)
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
core_x , core_y, core_z, L4b = DI.ConvertToBarycentricFrame(moon, core_pos, exportL4 = True, fixedL4 = True)

conf_x, conf_y, conf_z = co_x[0,:]-core_x[0], co_y[0,:]-core_y[0], co_z[0,:]-core_z[0]


print("core:")
print(np.round((core-L4[0])/1e3,4))
print("core_B:")
print(np.round((np.array([core_x[0], core_y[0], core_z[0]]).T -L4b[0])/1e3,4))
#%% 3D plot showing the orbit
RE = LL.Constants["Lunar_mean_dist"];
mu = LL.Constants["mu*"]
# L4b = np.ones(L4b.shape)*[0.5*RE, np.sqrt(3)/2 * RE, 0]

xlim, ylim, zlim = [-.4, 1.2], [-.1,1.3], [-.02, .02]
xticks = np.linspace(xlim[0],xlim[1],8)
yticks = np.linspace(ylim[0],ylim[1],8)
zticks = np.linspace(zlim[0],zlim[1],8)

fig,ax = PF.MakeFig(figsize = (10,10), dimensions='3D', proj = 'ortho')
# PF.Simple3DPlot(core_x[cof:],core_y[cof:],core_z[cof:],fig=fig,ax=ax, xmod=RE, ymod = RE, zmod = RE)
PF.Simple3DPlot(core_x[:cof],core_y[:cof],core_z[:cof],fig=fig,ax=ax, xmod=RE, ymod = RE, zmod = RE)

ax.scatter(-mu,0,0,c='r', marker='o', s=50)
ax.scatter(1-mu,0,0, c = 'b', marker = 'o')
ax.scatter(0.5, np.sqrt(3)/2,0,marker='d', c='k')

PF.DressFig(fig,ax, title= '', xlim= xlim, ylim = ylim, zlim = zlim,
            xlabel = 'x [R]', ylabel = 'y [R]', zlabel = 'z [R]',
            legendlist = ['Swarm motion - 1 year', 'Earth', 'Moon', 'L4' ])

ax.view_init(elev=90, azim=-90)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

ax.legend(['Swarm motion - 1 year', 'Earth', 'Moon', 'L4' ],loc=5)
fig.tight_layout()


ax.set_zticks([])
PF.saveFig(figure_folder,"fullyearorbit_z")

ax.view_init(elev=0, azim=-90)

ax.set_yticks([])
ax.set_zticks(zticks)
PF.saveFig(figure_folder,"fullyearorbit_y")

ax.view_init(elev=0, azim=0)

ax.set_xticks([])
ax.set_yticks(yticks)
PF.saveFig(figure_folder,"fullyearorbit_x")


fig,ax = PF.MakeFig(figsize = (10,10), dimensions='3D', proj = 'ortho')
# PF.Simple3DPlot(core_x[cof:],core_y[cof:],core_z[cof:],fig=fig,ax=ax, xmod=RE, ymod = RE, zmod = RE)
PF.Simple3DPlot(core_x,core_y,core_z,fig=fig,ax=ax, xmod=RE, ymod = RE, zmod = RE)

ax.scatter(-mu,0,0,c='r', marker='o', s=50)
ax.scatter(1-mu,0,0, c = 'b', marker = 'o')
ax.scatter(0.5, np.sqrt(3)/2,0,marker='d', c='k')

PF.DressFig(fig,ax, title= '', xlim= xlim, ylim = ylim, zlim = zlim,
            xlabel = 'x [R]', ylabel = 'y [R]', zlabel = 'z [R]',
            legendlist = ['Swarm motion - 5 years', 'Earth', 'Moon', 'L4' ])

ax.view_init(elev=90, azim=-90)

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_zticks(zticks)

fig.tight_layout()


ax.set_zticks([])
PF.saveFig(figure_folder,"fiveyearorbit_z")

ax.view_init(elev=0, azim=-90)

ax.set_yticks([])
ax.set_zticks(zticks)
PF.saveFig(figure_folder,"fiveyearorbit_y")

ax.view_init(elev=0, azim=0)

ax.set_xticks([])
ax.set_yticks(yticks)
PF.saveFig(figure_folder,"fiveyearorbit_x")

combinedOrbitPlot(core_x[:cof],core_y[:cof],core_z[:cof],"fullyearorbit_combined")

combinedOrbitPlot(core_x,core_y,core_z,"fiveyearorbit_combined",cof)
# PF.Simple3DPlot(core_pos[:,0],core_pos[:,1],core_pos[:,2], fig=fig,ax=ax)
# PF.Simple3DPlot(moon[:,0],moon[:,1],moon[:,2], fig=fig,ax=ax)

#%% Plot baselines up to the cutoff point
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x[:cof,:], BL_y[:cof,:], BL_z[:cof,:], includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m[:cof,:], 
                   baselineRateMagnitude = BLr_m[:cof,:], baselineRateMagnitudeThreshold = 1)
pp.savefig(figure_folder + "uvw_baselines")

#%%
coc = 2191*BL_x.shape[1]

# cstart = 3333 - 25
# cstop = 3333 + 25
# norm = np.array([[1,0,0]]).T
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[:coc],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img[cstart:cstop,cstart:cstop],
#                  '', True, figure_folder, 'PSF_x2',axes = 'psf')

# norm = np.array([[0,1,0]]).T
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[:coc],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img[cstart:cstop,cstart:cstop],
#                  '', True, figure_folder, 'PSF_y2',axes = 'psf')

# #%%
# norm = np.array([[0,0,1]]).T
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[:coc],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img[cstart:cstop,cstart:cstop],
#                  '', True, figure_folder, 'PSF_z2',axes = 'psf')
#%%

# #%%
# norm = np.array([[0,0,1]]).T
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[:coc],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img.T,
#                  '', True, figure_folder, 'PSF_z',axes = 'psf')

#%%
# c = 0
# orthogon = np.zeros((BL_m.shape[0],int(BL_m.shape[1]/2)))
# for i in range(0,Nsats-1):
#     for j in range(i+1,Nsats):
#         L = len(vx[:,i])
#         vi = np.concatenate([vx[:,i].reshape((L,1)),vy[:,i].reshape((L,1)), 
#                        vz[:,i].reshape((L,1))],axis=1)
#         vj = np.concatenate([vx[:,j].reshape((L,1)),vy[:,j].reshape((L,1)), 
#                        vz[:,j].reshape((L,1))],axis=1)
        
#         for k in range(0,L):
#             orthogon[k,c] = np.arccos(np.dot(vi[k],vj[k])/
#                                       (np.linalg.norm(vi[k])*np.linalg.norm(vj[k])))
            
#         c += 1
# #%%
# co = 2191
# td = (t-t[0])/(24*3600)
# fig,ax = PF.MakeFig(figsize = (12,3))

# for i in range(0,Nsats):
#     pp.plot(td[:co],orthogon[:co,i], alpha = 0.2, color='b')
# pp.plot(td[:co],np.max(orthogon[:co],axis=1), color = 'k')
# pp.plot(td[:co],np.min(orthogon[:co],axis=1), color = 'k')
# pp.axhline(np.pi/2,color='r')
# pp.yscale('Log')
# pp.xlabel('time in orbit [days]')
# pp.ylabel('$\omega$ [rad]')
# pp.xlim([0,365])
# pp.tight_layout()
# pp.savefig(figure_folder + f"angle_dot_product_{Nsats}")
#%% 
N_baseline_too_large = np.sum(BL_m[:cof,:] > 100e3)
N_baseline_too_small = np.sum(BL_m[:cof,:] < 500)
print( f"  this optima has {N_baseline_too_large} instances with too large baselines, \n \
      and {N_baseline_too_small} instances with too small baselines")
      
      
stab_lim_wanderoff = -1
stab_lim_minbaseline = -1
index = np.where(BL_m > 100e3)[0]
if len(index) != 0: stab_lim_wanderoff = np.min(index)
index = np.where(BL_m <500)[0]
if len(index) != 0: stab_lim_minbaseline = np.min(index)



print( f"  Upper stability limit of this orbit is {stab_lim_wanderoff*(4/24)} days untill a satellite wanders, \n \
      and {stab_lim_minbaseline*(4/24)} days untill a near-collision event")
      
    
# PF.UVWBaselinePlot(BL_x, BL_y, BL_z,
#                    includeNegatives=False, highlightfaulty = True, baselineMagnitude = BL_m, 
#              baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
# pp.savefig(figure_folder + 'uvw_baselines_1yr')
#%% baseline magnitude plot


fig,ax = PF.MakeFig(figsize=(12,4))
sz = int(BL_m.shape[1]/2)
yr = 365.25*24*3600
tl = t - t[0]
for j in range(0,sz):
    pp.plot(tl/yr,BL_m[:,j], color = 'b', alpha = 7.6/sz)
    
#legend = sz*[] + ['Near collision limit'] + ['Drift limit'] + ['Constellation maximum'] + ['Constellation minimum']
pp.axhline(500, color = 'r', ls='--', label='Near collision limit')
pp.axhline(1e3, color = 'r', ls='--', label='Soft detection limit')
pp.axhline(100e3, color = 'r',label='100 km drift limit')

pp.plot(tl/yr, np.max(BL_m,axis=1), color = 'k', label='Maximum baseline')
pp.plot(tl/yr, np.min(BL_m,axis=1), color = 'k', label='Minimum baseline')
pp.yscale('log')

pp.xlabel('time since 1 Jan 2030 [years]')
pp.ylabel("Baseline magnitude [m]")

pp.legend( loc=2)

pp.savefig(figure_folder + f"Baselinedistribution_{set_name}")


fig,ax = PF.MakeFig(figsize=(12,4))
sz = int(BL_m.shape[1]/2)
yr = 365.25*24*3600
tl = t - t[0]
for j in range(0,sz):
    pp.plot(tl/yr,BLr_m[:,j], color = 'r', alpha = 7.6/sz)
    
#legend = sz*[] + ['Near collision limit'] + ['Drift limit'] + ['Constellation maximum'] + ['Constellation minimum']
pp.axhline(1, color = 'r', label='1 m/s baseline rate')
# pp.axhline(100e3, color = 'r',label='100 km drift limit')

pp.plot(tl/yr, np.max(BLr_m,axis=1), color = 'k', label='Maximum baseline rate')
pp.plot(tl/yr, np.min(BLr_m,axis=1), color = 'k', label='Minimum baseline rate')
pp.yscale('log')

pp.xlabel('time since 1 Jan 2030 [years]')
pp.ylabel("Baseline rate [m/s]")

pp.legend( loc=2)

pp.savefig(figure_folder + f"Baselineratedistribution_{set_name}")

#%%
EM = 364522821.08369416
# fig,ax = PF.MakeFig(figsize=(8,8), dimensions = '3D')
# pp.plot(core_pos[:,0]/EM, core_pos[:,1]/EM, core_pos[:,2]/EM, label = 'Swarm position')
# ax.scatter(0,0,0,s = 30)
# pp.plot(moon[:,0]/EM, moon[:,1]/EM, moon[:,2]/EM, label = "Moon motion")

#%%
# test = sat[:,0] - core_pos[:,0]
# pp.figure()
# pp.plot(tl/yr,np.abs(test)/1000)
# pp.axhline(.5, color = 'r', ls='--', label='Near collision limit')
# pp.axhline(100, color = 'r',label='100 km drift limit')
# pp.yscale('log')

#%% Combined plot

fig,ax = PF.MakeFig(figsize=(16,12))
sz = int(BL_m.shape[1]/2)
yr = 365.25*24*3600
tl = t - t[0]
for j in range(0,sz):
    pp.plot(tl/yr,BL_m[:,j], color = 'b', alpha = 7.6/sz)
    pp.plot(tl/yr,BLr_m[:,j], color = 'r', alpha = 7.6/sz)
    
#legend = sz*[] + ['Near collision limit'] + ['Drift limit'] + ['Constellation maximum'] + ['Constellation minimum']
pp.axhline(500, color = 'r', ls='--', label='Near collision limit')
pp.axhline(100e3, color = 'r',label='100 km drift limit')

pp.plot(tl/yr, np.max(BL_m,axis=1), color = 'k', label='Maximum baseline')
pp.plot(tl/yr, np.min(BL_m,axis=1), color = 'k', label='Minimum baseline')
pp.plot(tl/yr, np.mean(BL_m,axis=1), color = 'k',ls='--', label='Mean baseline')

pp.plot(tl/yr, np.max(BLr_m,axis=1), color = 'g', label='Maximum baseline rate')
pp.plot(tl/yr, np.min(BLr_m,axis=1), color = 'g', label='Minimum baseline rate')
pp.plot(tl/yr, np.mean(BLr_m,axis=1), color = 'g',ls='--', label='Mean baseline rate')
pp.axhline(1, color = 'r', label='1 m/s baseline rate')
pp.yscale('log')

pp.xlabel('time since 1 Jan 2030 [years]')
pp.ylabel("Baseline magnitude [m] / rate [m/s]")

pp.legend( loc=2)
pp.tight_layout()
pp.savefig(figure_folder + f"combinedplot_{set_name}")

#%%
co = 2400
yr = 24*3600
fig,ax = PF.MakeFig(figsize=(12,4))
for j in range(0,sz):
    pp.plot(tl[:co]/yr,BL_m[:co,j], color = 'b', alpha = 7.6/sz)
    pp.plot(tl[:co]/yr,BLr_m[:co,j], color = 'g', alpha = 7.6/sz)
    
#legend = sz*[] + ['Near collision limit'] + ['Drift limit'] + ['Constellation maximum'] + ['Constellation minimum']
pp.axhline(500, color = 'r', ls='--', label='Near collision limit')
pp.axhline(100e3, color = 'r',label='100 km drift limit')

pp.plot(tl[:co]/yr, np.max(BL_m[:co],axis=1), color = 'k', label='Maximum baseline')
pp.plot(tl[:co]/yr, np.min(BL_m[:co],axis=1), color = 'k', label='Minimum baseline')
pp.plot(tl[:co]/yr, np.mean(BL_m[:co],axis=1), color = 'y', label='Mean baseline')

pp.plot(tl[:co]/yr, np.max(BLr_m[:co],axis=1), color = 'm', label='Maximum baseline rate')
pp.plot(tl[:co]/yr, np.min(BLr_m[:co],axis=1), color = 'm', label='Minimum baseline rate')
pp.plot(tl[:co]/yr, np.mean(BLr_m[:co],axis=1), color = 'm',ls='--', label='Mean baseline rate')
pp.axhline(1, color = 'r', label='1 m/s baseline rate')
pp.yscale('log')

pp.xlabel('time since 1 Jan 2030 [days]')
pp.ylabel("Baseline magnitude [m] / rate [m/s]")

pp.legend( loc=2, bbox_to_anchor=(1.05, 1))
pp.xlim([0, tl[co]/yr])
pp.tight_layout()
pp.savefig(figure_folder + f"combinedplot_1yr_{set_name}")


co = len(tl)
yr *=365.25
fig,ax = PF.MakeFig(figsize=(12,4))
for j in range(0,sz):
    pp.plot(tl[:co]/yr,BL_m[:co,j], color = 'b', alpha = 7.6/sz)
    pp.plot(tl[:co]/yr,BLr_m[:co,j], color = 'g', alpha = 7.6/sz)
    
#legend = sz*[] + ['Near collision limit'] + ['Drift limit'] + ['Constellation maximum'] + ['Constellation minimum']
pp.axhline(500, color = 'r', ls='--', label='Near collision limit')
pp.axhline(100e3, color = 'r',label='100 km drift limit')

pp.plot(tl[:co]/yr, np.max(BL_m[:co],axis=1), color = 'k', label='Maximum baseline')
pp.plot(tl[:co]/yr, np.min(BL_m[:co],axis=1), color = 'k', label='Minimum baseline')
pp.plot(tl[:co]/yr, np.mean(BL_m[:co],axis=1), color = 'y', label='Mean baseline')

pp.plot(tl[:co]/yr, np.max(BLr_m[:co],axis=1), color = 'm', label='Maximum baseline rate')
pp.plot(tl[:co]/yr, np.min(BLr_m[:co],axis=1), color = 'm', label='Minimum baseline rate')
pp.plot(tl[:co]/yr, np.mean(BLr_m[:co],axis=1), color = 'm',ls='--', label='Mean baseline rate')
pp.axhline(1, color = 'r', label='1 m/s baseline rate')
pp.yscale('log')

pp.xlabel('time since 1 Jan 2030 [years]')
pp.ylabel("Baseline magnitude [m] / rate [m/s]")

pp.legend( loc=2, bbox_to_anchor=(1.05, 1))
pp.xlim([0, tl[-1]/yr])
pp.tight_layout()
pp.savefig(figure_folder + f"combinedplot_5yr_{set_name}")


#%%
ang = np.cross(moon[0,:3],moon[0,3:])/(np.linalg.norm(moon[0,:3])**2)
VL4 = np.cross(ang,L4[0])

rel_vx_0 = vx[0] - VL4[0]
rel_vy_0 = vy[0] - VL4[1]
rel_vz_0 = vz[0] - VL4[2]

H  = np.cross(moon[0,:3],moon[0,3:])
b_x = moon[0,:3]/np.linalg.norm(moon[0,:3])
b_z = H/np.linalg.norm(H)
b_y = np.cross(b_z,b_x)

rel_v = np.array([rel_vx_0[0], rel_vy_0[0], rel_vz_0[0]])

b_velx = vectorProj(rel_v,b_x)
b_vely = vectorProj(rel_v,b_y)
b_velz = vectorProj(rel_v,b_z)

b_vel = np.array([b_velx, b_vely, b_velz])

fig,ax = PF.MakeFig(dimensions = '3D')
pp.title('barycentric frame swarm relative to core')
# pp.plot(core_x[:2].flatten()-L4b[:2,0],core_y[:2].flatten()-L4b[:2,1],
#         core_z[:2].flatten()-L4b[:2,2])
ax.scatter((co_x[0]-L4b[0,0])/1e3, (co_y[0]-L4b[0,1])/1e3, (co_z[0]-L4b[0,2])/1e3)
ax.scatter(0,0,0,s = 60, marker = 'd', color = 'g')

ax.scatter((core_x[0]-L4b[0,0])/1e3, (core_y[0]-L4b[0,1])/1e3, (core_z[0]-L4b[0,2])/1e3,marker='s',color='r')
ax.quiver((core_x[0]-L4b[0,0])/1e3, (core_y[0]-L4b[0,1])/1e3, (core_z[0]-L4b[0,2])/1e3, b_vel[0],
          b_vel[1], b_vel[2],length=1.2)

ax.set_xlabel('$x_B$ [km]')
ax.set_ylabel('$y_B$ [km]')
ax.set_zlabel('$z_B$ [km]')

pp.savefig(figure_folder + '3d-orbit-barycentric')
#%%

# xlim, ylim, zlim = [-.4, 1.2], [-.1,1.3], [-2, 2]
    
fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.scatter((co_x[0]-L4b[0,0])/1e3, (co_y[0]-L4b[0,1])/1e3,)
ax1.scatter(0,0,s = 60, marker = 'd', color = 'g')
ax1.scatter((core_x[0]-L4b[0,0])/1e3, (core_y[0]-L4b[0,1])/1e3,marker='s',color='r')
ax1.quiver((core_x[0]-L4b[0,0])/1e3, (core_y[0]-L4b[0,1])/1e3, b_vel[0],
          b_vel[1])

# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')

ax3.scatter( (co_y[0]-L4b[0,1])/1e3, (co_z[0]-L4b[0,2])/1e3)
ax3.scatter(0,0,s = 60, marker = 'd', color = 'g')
ax3.scatter( (core_y[0]-L4b[0,1])/1e3, (core_z[0]-L4b[0,2])/1e3,marker='s',color='r')
ax3.quiver( (core_y[0]-L4b[0,1])/1e3, (core_z[0]-L4b[0,2])/1e3, 
          b_vel[1], b_vel[2])

# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')

ax2.scatter((co_x[0]-L4b[0,0])/1e3,  (co_z[0]-L4b[0,2])/1e3)
ax2.scatter(0,0,s = 60, marker = 'd', color = 'g')
ax2.scatter((core_x[0]-L4b[0,0])/1e3, (core_z[0]-L4b[0,2])/1e3,marker='s',color='r')
ax2.quiver((core_x[0]-L4b[0,0])/1e3,  (core_z[0]-L4b[0,2])/1e3, b_vel[0],
          b_vel[2])
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')

fig.legend(['Satellites', 'L4', 'Core position', 'Core velocity' ],bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_visualisation")

#%%

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.scatter((co_x[0]-core_x[0])/1e3, (co_y[0]-core_y[0])/1e3,)
ax1.scatter(0, 0,marker='s',color='r')
ax1.quiver(np.zeros(b_vel[1].shape), np.zeros(b_vel[1].shape), b_vel[0],
          b_vel[1])

# ax1.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_y[:10].flatten()-L4b[:10,1])/1e3)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')
ax1.set_xlim([-50,50])
ax1.set_ylim([-50,50])

ax3.scatter( (co_y[0]-core_y[0])/1e3, (co_z[0]-core_z[0])/1e3)
ax3.scatter( 0,0,marker='s',color='r')
ax3.quiver( np.zeros(b_vel[1].shape),np.zeros(b_vel[1].shape), 
          b_vel[1], b_vel[2])

# ax3.plot( (core_y[:10].flatten()-L4b[:10,1])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')
ax3.set_xlim([-50,50])
ax3.set_ylim([-50,50])

ax2.scatter((co_x[0]-core_x[0])/1e3,  (co_z[0]-core_z[0])/1e3)
ax2.scatter(0,0,marker='s',color='r')
ax2.quiver(np.zeros(b_vel[1].shape),np.zeros(b_vel[1].shape), b_vel[0],
          b_vel[2])
# ax2.plot( (core_x[:10].flatten()-L4b[:10,0])/1e3,(core_z[:10].flatten()-L4b[:10,2])/1e3)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')
ax2.set_xlim([-50,50])
ax2.set_ylim([-50,50])

fig.legend(['Satellites', 'Core position', 'Core velocity' ],bbox_to_anchor=(0.65, 0.95))
pp.tight_layout()
pp.savefig(figure_folder + "swarm_visualisation_core")


#%%

# fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

# ax1.scatter((co_x[0]-L4b[0,0])/1e3, (co_y[0]-L4b[0,1])/1e3,)
# ax1.scatter(0,0,s = 60, marker = 'd', color = 'g')
# ax1.scatter((core_x[0]-L4b[0,0])/1e3, (core_y[0]-L4b[0,1])/1e3,marker='s',color='r')
# ax1.quiver((core_x[0]-L4b[0,0])/1e3, (core_y[0]-L4b[0,1])/1e3, b_vel[0],
#           b_vel[1])

# ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax1.set_xlabel('$x_B$ [km]')
# ax1.set_ylabel('$y_B$ [km]')

# ax3.scatter( (co_y[0]-L4b[0,1])/1e3, (co_z[0]-L4b[0,2])/1e3)
# ax3.scatter(0,0,s = 60, marker = 'd', color = 'g')
# ax3.scatter( (core_y[0]-L4b[0,1])/1e3, (core_z[0]-L4b[0,2])/1e3,marker='s',color='r')
# ax3.quiver( (core_y[0]-L4b[0,1])/1e3, (core_z[0]-L4b[0,2])/1e3, 
#           b_vel[1], b_vel[2])


# ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax3.set_xlabel('$y_B$ [km]')
# ax3.set_ylabel('$z_B$ [km]')

# ax2.scatter((co_x[0]-L4b[0,0])/1e3,  (co_z[0]-L4b[0,2])/1e3)
# ax2.scatter(0,0,s = 60, marker = 'd', color = 'g')
# ax2.scatter((core_x[0]-L4b[0,0])/1e3, (core_z[0]-L4b[0,2])/1e3,marker='s',color='r')
# ax2.quiver((core_x[0]-L4b[0,0])/1e3,  (core_z[0]-L4b[0,2])/1e3, b_vel[0],
#           b_vel[2])

# ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax2.set_xlabel('$x_B$ [km]')
# ax2.set_ylabel('$z_B$ [km]')

# for i in range(0,Nsats):
    
#     ax1.plot( (co_x[:10,i]-L4b[:10,0])/1e3,(co_y[:10,i]-L4b[:10,1])/1e3,color='r', alpha = 0.3)
#     ax2.plot( (co_x[:10,i]-L4b[:10,0])/1e3,(co_z[:10,i]-L4b[:10,2])/1e3,color='r', alpha = 0.3)
#     ax3.plot( (co_y[:10,i]-L4b[:10,1])/1e3,(co_z[:10,i]-L4b[:10,2])/1e3,color='r', alpha = 0.3)

# fig.legend(['Satellites', 'L4', 'Core position', 'Core velocity' ],bbox_to_anchor=(0.65, 0.95))
# pp.tight_layout()
# pp.savefig(figure_folder + "swarm_visualisation_with_orbits")

#%%
co = int(365.25*24/4)
xlim, ylim, zlim = [-110, 110], [-110,110], [-110, 110]

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))


ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('x [km]')
ax1.set_ylabel('y [km]')
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)

ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('y [km]')
ax3.set_ylabel('z [km]')
ax3.set_xlim(ylim)
ax3.set_ylim(zlim)

ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('x [km]')
ax2.set_ylabel('z [km]')
ax2.set_xlim(xlim)
ax2.set_ylim(zlim)

for i in range(0,2*sz):
    ax1.plot(BL_x[:co,i]/1e3,BL_y[:co,i]/1e3,color='b')
    ax1.plot(BL_x[:co,i]/1e3,BL_y[:co,i]/1e3,color='r', alpha=5/sz,zorder = 10)
    ax2.plot(BL_y[:co,i]/1e3,BL_z[:co,i]/1e3,color='b')
    ax2.plot(BL_y[:co,i]/1e3,BL_z[:co,i]/1e3,color='r', alpha=5/sz,zorder = 10)
    ax3.plot(BL_x[:co,i]/1e3,BL_z[:co,i]/1e3,color='b')
    ax3.plot(BL_x[:co,i]/1e3,BL_z[:co,i]/1e3,color='r', alpha=5/sz,zorder = 10)

pp.tight_layout()
pp.savefig(figure_folder + "baseline_profile_1yr")


#%%

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))

ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('x [km]')
ax1.set_ylabel('y [km]')
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)

ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('y [km]')
ax3.set_ylabel('z [km]')
ax3.set_xlim(ylim)
ax3.set_ylim(zlim)

ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('x [km]')
ax2.set_ylabel('z [km]')
ax2.set_xlim(xlim)
ax2.set_ylim(zlim)
sel = BL_m <100e3
for i in range(0,2*sz):
    ax1.plot(BL_x[sel[:,i],i]/1e3,BL_y[sel[:,i],i]/1e3,color='b')
    ax1.plot(BL_x[sel[:,i],i]/1e3,BL_y[sel[:,i],i]/1e3,color='r', alpha=5/sz,zorder = 10)
    ax2.plot(BL_y[sel[:,i],i]/1e3,BL_z[sel[:,i],i]/1e3,color='b')
    ax2.plot(BL_y[sel[:,i],i]/1e3,BL_z[sel[:,i],i]/1e3,color='r', alpha=5/sz,zorder = 10)
    ax3.plot(BL_x[sel[:,i],i]/1e3,BL_z[sel[:,i],i]/1e3,color='b')
    ax3.plot(BL_x[sel[:,i],i]/1e3,BL_z[sel[:,i],i]/1e3,color='r', alpha=5/sz,zorder = 10)

pp.tight_layout()
pp.savefig(figure_folder + "baseline_profile_5yr")

# #%%
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[sel.flatten()],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img[cstart:cstop,cstart:cstop],
#                   '', True, figure_folder, 'PSF_x_ULT',axes = 'psf')

# norm = np.array([[0,1,0]]).T
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[sel.flatten()],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img[cstart:cstop,cstart:cstop],
#                   '', True, figure_folder, 'PSF_y_ULT',axes = 'psf')

# #%%
# norm = np.array([[0,0,1]]).T
# I_PSF, I_S,cost = GF.PSFanalysis(norm, BL_c[sel.flatten()],resolution = 30)
# img = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,0])), np.uint8)
# PF.Simple2DImage(img[cstart:cstop,cstart:cstop],
#                   '', True, figure_folder, 'PSF_z_ULT',axes = 'psf')
#%%
fig,ax = PF.MakeFig(dimensions = '3D')
ax.quiver(0,0,0,b_x[0],b_x[1],b_x[2])
ax.quiver(0,0,0,b_y[0],b_y[1],b_y[2])
ax.quiver(0,0,0,b_z[0],b_z[1],b_z[2])

ax.scatter(0,0,0)
m = moon[0,:3]/np.linalg.norm(moon[0,:3])
ax.scatter(m[0],m[1],m[2])
l = L4/np.linalg.norm(moon[0,:3])

ax.scatter(l[0,0],l[0,1],l[0,2])

#%% print data

xcc = np.round((x[0]-core[0])/1e3,4)
ycc = np.round((y[0]-core[1])/1e3,4)
zcc = np.round((z[0]-core[2])/1e3,4)
    
print("core velocity:")
print(np.round(rel_v,4))
print("x:")
print(xcc)
print("y:")
print(ycc)
print("z:")
print(zcc)



for i in range(0,Nsats):

    print(f"{xcc[i]} & {ycc[i]} & {zcc[i]}")
    
#%% soft detection limit
    
dt = 4*3600

vrel = 0.10

L = np.sqrt(2*(0.5*vrel*dt)**2)

print(f"system soft detection is {L} meters for a threshold of {vrel} m/s")

#%% plot blr versus blm
thres = 6.944e-2
fig,ax = PF.MakeFig(figsize=(12,3))

pp.scatter(BLr_m[BL_m < 1e3].flatten(),BL_m[BL_m < 1e3].flatten(), alpha = 0.4)
pp.scatter(BLr_m[BL_m < 1e3].flatten(),BL_m[BL_m < 1e3].flatten(), alpha = 0.01, color = 'r',zorder=10)
pp.xscale('log')    
pp.yscale('log')
pp.xlabel("baseline rate")
pp.ylabel("baseline magnitude")
pp.xlim([1e-3,10])
pp.ylim([200,1100])

pp.axvline(thres,color = 'r')
pp.tight_layout()
pp.savefig(figure_folder + f"baselinerate_vs_magnitude_{Nsats}_sat_focused")
sel = BLr_m[BL_m < 1e3].flatten()
L = len(sel)

lfail = len(sel[sel > thres].flatten())

print(f"Of {L} smaller than 1km baseline occurances, {lfail} has a velocity larger than {thres} cm/s")

print(f"percentage wise: {lfail/L * 100} %")

#%%

blm = BL_m[:co]
blr = BLr_m[:co]
fig,ax = PF.MakeFig(figsize=(12,3))
pp.scatter(blr[blm < 1e3].flatten(),blm[blm < 1e3].flatten(), alpha = 0.4)
pp.scatter(blr[blm < 1e3].flatten(),blm[blm < 1e3].flatten(), alpha = 0.01, color = 'r',zorder=10)
pp.xscale('log')    
pp.yscale('log')
pp.xlabel("baseline rate")  
pp.ylabel("baseline magnitude")
pp.xlim([1e-3,5])
pp.ylim([200,1100])

pp.axvline(thres,color = 'r')
pp.tight_layout()
pp.savefig(figure_folder + f"baselinerate_vs_magnitude_{Nsats}_sat_focused_1yr")
sel = blr[blm < 1e3].flatten()
L = len(sel)

lfail = len(sel[sel > thres].flatten())

print(f"Of {L} smaller than 1km baseline occurances, {lfail} has a velocity larger than {thres} cm/s")

print(f"percentage wise: {lfail/L * 100} %")

#%% Barycentric baseline plots

BL_x, BL_y, BL_z, BL_c, BL_m = GF.GetBaseLinesNoV(co_x,co_y,co_z, t)

co = int(365.25*24/4)
xlim, ylim, zlim = [-110, 110], [-110,110], [-110, 110]

fig, (ax1,ax2,ax3) = pp.subplots(1,3, figsize = (15,5))


ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.set_xlabel('$x_B$ [km]')
ax1.set_ylabel('$y_B$ [km]')
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)

ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.set_xlabel('$y_B$ [km]')
ax3.set_ylabel('$z_B$ [km]')
ax3.set_xlim(ylim)
ax3.set_ylim(zlim)

ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.set_xlabel('$x_B$ [km]')
ax2.set_ylabel('$z_B$ [km]')
ax2.set_xlim(xlim)
ax2.set_ylim(zlim)

for i in range(0,2*sz):
    ax1.plot(BL_x[:co,i]/1e3,BL_y[:co,i]/1e3,color='b')
    ax1.plot(BL_x[:co,i]/1e3,BL_y[:co,i]/1e3,color='r', alpha=5/sz,zorder = 10)
    ax2.plot(BL_y[:co,i]/1e3,BL_z[:co,i]/1e3,color='b')
    ax2.plot(BL_y[:co,i]/1e3,BL_z[:co,i]/1e3,color='r', alpha=5/sz,zorder = 10)
    ax3.plot(BL_x[:co,i]/1e3,BL_z[:co,i]/1e3,color='b')
    ax3.plot(BL_x[:co,i]/1e3,BL_z[:co,i]/1e3,color='r', alpha=5/sz,zorder = 10)

pp.tight_layout()
pp.savefig(figure_folder + "baseline_profile_1yr_B")
#%% swarm fold plot
td = (t-t[0])/(24*3600)

conf_x, conf_y, conf_z = co_x-core_x, co_y-core_y, co_z-core_z
fig, (ax1,ax2,ax3) = pp.subplots(3,1, figsize = (15,7))

for i in range(0,15):
    ax1.plot(td[:co],conf_x[:co,i]/1e3,color='b', alpha = 0.5)
    ax2.plot(td[:co],conf_y[:co,i]/1e3,color='b', alpha = 0.5)
    ax3.plot(td[:co],conf_z[:co,i]/1e3,color='b', alpha = 0.5)

ax1.set_ylabel('Relative $x_B$[km]')
ax2.set_ylabel('Relative $y_B$[km]')

ax3.set_xlabel('time in orbit [days]')
ax3.set_ylabel('Relative $z_B$[km]')

ax1.set_xlim([0,365])
ax2.set_xlim([0,365])
ax3.set_xlim([0,365])

pp.tight_layout()
pp.savefig(figure_folder + 'harmonic_decay')