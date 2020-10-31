# -*- coding: utf-8 -*-
"""
Testing ground to try out the state vector propagation
Created on Tue Apr 21 13:12:38 2020

@author: Jurriez
"""
import numpy as np
import sys
from os import chdir
from os import getcwd

from scipy.integrate import odeint

sys.path.append('../')
chdir('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF
from matplotlib import pyplot as pp
chdir('Baseline Propagation')

def propagateState(state, dt):
    # More readable approach
    final = np.zeros(state.shape)
    x,y,z,dx,dy,dz,ddx,ddy,ddz = state
    
    # update position
    final[0] = x + dt*dx + 0.5*ddx*dt**2
    final[1] = y + dt*dy + 0.5*ddy*dt**2
    final[2] = z + dt*dz + 0.5*ddz*dt**2
    # update velocities
    final[3] = dx + dt*ddx
    final[4] = dy + dt*ddy
    final[5] = dz + dt*ddz
    # update acceleration with current state info
    final[6] =  2*state[4]
    final[7] = - 2*state[3]
    final[8] = - state[2]
    
    return final

def propagateFull(state,dt,mu):
    # Create copy arrays to avoid editing the same vector during integration
    orig = np.copy(state)
    final = np.zeros(state.shape)
    # Bind shortcuts to make the code more readable
    x,y,z,dx,dy,dz,ddx,ddy,ddz = orig
    
    
    # turn mu into mu*
    mu = 1-mu
    
    # Compute range vectors
    r1 = np.sqrt((x+mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x-1 +mu)**2 + y**2 + z**2)
    
    
    # update position
    final[0] = x + dt*dx + 0.5*ddx*dt**2 #x
    final[1] = y + dt*dy + 0.5*ddy*dt**2 #y
    final[2] = z + dt*dz + 0.5*ddz*dt**2 #z
    
    # update velocities
    final[3] = (dx + dt*ddx)
    final[4] = (dy + dt*ddy)
    final[5] = (dz + dt*ddz)

    # update acceleration with current state info
    final[6] = (- (1-mu)*(final[0] + mu )/r1**3 - mu*(final[0]-1+mu)/r2**3 +2*final[4] + final[0])
    final[7] = (- (1-mu)*final[1]/r1**3 -mu*final[1]/r2**3 - 2*final[3] + final[1])
    final[8] = (- (1 - mu)*final[2]/r1**3 - mu*final[2]/r2**3)
    
    return final

def propagateFullRK4(state,dt,mu):
    orig = np.copy(state)
    
    k1 = propagateRK4(orig,0,mu)
    k2 = propagateRK4(orig+dt/2*k1,dt/2,mu)
    k3 = propagateRK4(orig+dt/2*k2,dt/2,mu)
    k4 = propagateRK4(orig+dt/2*k3,dt,mu)
    
    final = orig + dt/6*(k1+2*k2+2*k3+k4)
    final = propagateFull(final,0,mu)
    
    
    return final

def propagateRK4(state,dt,mu):
    # Create copy arrays to avoid editing the same vector during integration
    orig = np.copy(state)
    final = np.zeros(state.shape)
    # Bind shortcuts to make the code more readable
    x,y,z,dx,dy,dz,ddx,ddy,ddz = orig
    
    
    # turn mu into mu*
    mu = 1-mu
    
    # Compute range vectors
    r1 = np.sqrt((x+mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x-1 +mu)**2 + y**2 + z**2)
    
    
    # update position
    final[0] = dx
    final[1] = dy
    final[2] = dz
    
    # update velocities
    final[3] = ddx
    final[4] = ddy
    final[5] = ddz

    # update acceleration with current state info
    final[6] = (- (1-mu)*(final[0] + mu )/r1**3 - mu*(final[0]-1+mu)/r2**3 +2*final[4] + final[0]) - ddx
    final[7] = (- (1-mu)*final[1]/r1**3 -mu*final[1]/r2**3 - 2*final[3] + final[1]) - ddy
    final[8] = (- (1 - mu)*final[2]/r1**3 - mu*final[2]/r2**3) - ddz
    
    return final

def psuedopot(x,y,mu,r1,r2):
    return (1-mu)/r1 + mu/r2 + (x**2 + y**2)/2
    

G = 6.674e-11
m_E = 5.97237e24
m_M = 7.34767309e22
m_t = m_E + m_M
mu = m_E/(m_t)


L4 = np.zeros(9)
L4[0] = mu-0.5
L4[1] = np.sqrt(3)/2

#%% Import tudat reference data 

chdir("../")
cwd = getcwd()
append = cwd + "/Data/"
print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")
#%%

datafolder = cwd + '/Data/PerturbationAnalysis/'

addendum = "_365d"
filename = datafolder + "propagationHistory_0_30km.dat"
t2, xd, yd, zd, vxd, vyd, vzd, md  = DI.ImportPropagationHistory(filename,1, True)

EM = np.sqrt(md[0,0]**2 + md[0,1]**2 + md[0,2]**2) #364522821.08369416 # Tudat t_0

state_s1 = np.concatenate((xd[:,:1],yd[:,:1],zd[:,:1],vxd[:,:1],vyd[:,:1],vzd[:,:1]),axis=1)
pos_s1, vel_s1 = DI.convertFullStateToBarycentric(md,state_s1)
pos_s1 /= EM
vel_s1 /= EM

jjk1 = pos_s1 - L4[:3]

state_s2 = np.concatenate((xd[:,1:],yd[:,1:],zd[:,1:],vxd[:,1:],vyd[:,1:],vzd[:,1:]),axis=1)
pos_s2, vel_s2 = DI.convertFullStateToBarycentric(md,state_s2)
pos_s2 /= EM
vel_s2 /= EM

jjk = pos_s2 - L4[:3]

#%%
char_time = np.sqrt(EM**3/(G*m_t))
vel_s1 *= char_time
vel_s2 *= char_time


S1,S2 = ( np.copy(L4) for i in range(0,2))
S1[0:3] = pos_s1[0,:]
S1[3:6] = vel_s1[0,:]
S2[0:3] = pos_s2[0,:]
S2[3:6] = vel_s2[0,:]
# S2[1] += 1e2/EM
# S2[0] += 100 /EM
# S2[1] += 0 /EM
# # S2[2] += 0 /EM
# S2[3] += 1 *char_time/EM # need to make m/s dimensionless, thus *tc/EM
# S2[4] += 0 *char_time/EM
# S2[5] += 11 *char_time/EM

S1f = np.copy(S1)
S2f = np.copy(S2)

S1 = propagateState(S1,0)
S2 = propagateState(S2,0)
S1f = propagateFull(S1f,0,mu)
S2f = propagateFull(S2f,0,mu)
S2f0 = np.copy(S2f)
dS = S1-S2
dSf = S1f-S2f

t0 = 0
t1 = 365.25*24*3600
dt = 30*60
N = np.int((t1-t0)/dt)+1
t = np.zeros(N)


c_S1, c_S2, c_dS,c_L1,c_L2 = ( np.zeros((N,9)) for i in range(0,5))
c_S1f, c_S2f, c_dSf,c_L1f,c_L2f = ( np.zeros((N,9)) for i in range(0,5))

#%% Run a manual solver

# Correct timestep to account for characteristic time of system
dt /=char_time
for i in range(0,N):
    t[i] = i*dt
    c_S1[i,:],c_S2[i,:],c_dS[i,:] = S1,S2,dS
    c_L1[i,:],c_L2[i,:] = S1-L4,S2-L4
    
    c_S1f[i,:],c_S2f[i,:],c_dSf[i,:] = S1f,S2f,dSf
    c_L1f[i,:],c_L2f[i,:] = S1f-L4,S2f-L4
    
    S1 = propagateState(S1,dt)
    S2 = propagateState(S2,dt)
    S1f = propagateFullRK4(S1f,dt,mu)
    S2f = propagateFullRK4(S2f,dt,mu)
    dS = S1-S2
    dSf = S1-S2
    
error_fullmodel = pos_s1 - c_S1f[:,:3]

#%%
legend = ['x','y', 'z', 'dx', 'dy', 'dz', 'ddx', 'ddy', 'ddz']
PF.ClearPlots()
fig,ax = PF.Matrix2DPlot(t,np.abs(c_S2),xmod =1 , ymod= 1, xlabel = 't[s]',
                         logScale=True, legendlist = legend)

#%% 

scale = 1e3/EM

fig,ax = PF.Simple3DPlot(c_L2[:,0],c_L2[:,1],c_L2[:,2],
                         xmod=scale,ymod=scale, zmod=scale)
PF.Simple3DScatter(0,0,0,title='Euler propagated orbit', fig=fig,ax=ax,
                   marker='*',color='r')
# ax.set_xlim([-0,0.1])
# ax.set_ylim([-0.1,0.1])
# ax.set_zlim([-0.1,0.1])
PF.saveFig(LL.folder_dict["NumericalSim"], "3bp orbit L4")

fig,ax = PF.Simple3DPlot(c_L2f[:,0],c_L2f[:,1],c_L2f[:,2],
                         xmod=scale,ymod=scale, zmod=scale)
PF.Simple3DScatter(0,0,0,title='RK4 propagated orbit', 
                   fig=fig,ax=ax, marker='*',color='r')
# ax.set_zlim([-0.1,0.1])
PF.saveFig(LL.folder_dict["NumericalSim"], "RK4_prop_orbit")

fig,ax = PF.Simple3DPlot(c_S2f[:,0],c_S2f[:,1],c_S2f[:,2],
                         xmod=scale,ymod=scale, zmod=scale)
PF.Simple3DScatter(0,0,0,title='Satellite S2 full model orbit relative to Earth', 
                   fig=fig,ax=ax, marker='*',color='r')
# ax.set_zlim([-1,1])


#%% Scatter plot with planet positions and L4 locations
fig,ax = PF.Simple3DPlot(jjk[:,0],jjk[:,1],
                         jjk[:,2], 
                         xmod=(1000/EM),ymod=1000/EM, zmod=1000/EM)
# PF.Simple3DScatter(0,0,0,title='Full perturbation model orbit around L4',
#                    fig=fig,ax=ax, marker='*',color='r')
pp.title('Full perturbation model orbit around L4')
ax.scatter(0,0,0,s = 200, marker = '*', color = 'r')
PF.saveFig(LL.folder_dict["NumericalSim"], "full model orbit around L4")

BL = pos_s1 - pos_s2

PF.Simple3DPlot(BL[:,0], BL[:,1], BL[:,2], title='Baseline motion during orbit - full perturbation model',
                savefig = False, figfolder = '',
                name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]',
                zlabel='z [km]', xmod = scale, ymod = scale, zmod = scale,
                xlim = None, ylim = None , zlim = None, fig = None, ax = None)

# l4 = DI.computeL4Location(md,mu)
# PF.PlotAroundEarth(state_s1[:,0].reshape((N,1)),state_s1[:,1].reshape((N,1))
#                    ,state_s1[:,2].reshape((N,1)), md, l4, plotMoonTrails = False)
#%%
fig,ax = PF.Simple3DScatter(pos_s2[0,0],pos_s2[0,1],pos_s2[0,2], xmod=1,ymod=1,zmod=1)
PF.Simple3DScatter(c_S2[0,0],c_S2[0,1],c_S2[0,2], fig=fig,ax=ax, xmod=1,ymod=1,zmod=1,marker='*')
PF.Simple3DScatter(L4[0],L4[1],L4[2], fig=fig,ax=ax, xmod=1,ymod=1,zmod=1)
PF.Simple3DScatter(0,0,0, fig=fig,ax=ax, xmod=1,ymod=1,zmod=1,marker='d')
PF.Simple3DScatter(1,0,0, fig=fig,ax=ax, xmod=1,ymod=1,zmod=1,marker='d')
# ax.set_zlim([-1,1])

#%% Initials state showcase

fig,ax = PF.MakeFig(figsize=(12,12), dimensions='3D')
PF.Simple3DScatter(0,0,0,xmod = 1, ymod=1, zmod = 1, color='b',fig=fig,ax=ax)
PF.Simple3DScatter(md[0,0]/EM,md[0,1]/EM,md[0,2]/EM,xmod = 1, ymod=1, zmod = 1,
                   color='b',fig=fig,ax=ax)

p0 = md[0,:3]/EM 
p1 = p0 + 0.1*md[0,3:]/np.linalg.norm(md[0,3:])
PF.AddVector(ax, p0, p1)

PF.Simple3DScatter(state_s1[0,0]/EM,state_s1[0,1]/EM,state_s1[0,2]/EM,xmod = 1, ymod=1, zmod = 1,
                   color='k',fig=fig,ax=ax,marker='d')
PF.Simple3DScatter(state_s2[0,0]/EM,state_s2[0,1]/EM,state_s2[0,2]/EM,xmod = 1, ymod=1, zmod = 1,
                   color='k',fig=fig,ax=ax,marker='d')

#%% Compare full model orbit to reference orbit
# initial velocities are wrong!
error_3bp = c_S1f[:,:3] - pos_s1
error_3bp *=EM

relpos_orig = pos_s1 - pos_s2
relpos_prop = c_S1f[:,:3] - c_S2f[:,:3]

error_baseline_fullmodel = (relpos_prop - relpos_orig)*EM

tp = t*char_time/(24*3600)
N = int(50*24*3600/1800)
fig,ax = PF.MakeFig(figsize=(8,8))

PF.Simple2DPlot(tp[:N], np.linalg.norm(error_3bp[:N,:],axis=1),
                xmod=1, ymod = 1, fig = fig, ax = ax)
PF.Simple2DPlot(tp[:N], np.linalg.norm(error_baseline_fullmodel[:N,:],axis=1),
                title='Analytical 3 body problem baseline error compared to full perturbation model',
                xlabel='Time [days]',
                ylabel = 'Baseline error [m]', xmod=1, ymod = 1, fig = fig, ax = ax)

ax.set_yscale("Log")
ax.set_ylim([1e-1, 1e8])
ax.set_xticks(np.arange(min(tp), tp[N]+10, 10))
ax.legend(["Position error", "Baseline error"])

PF.saveFig(LL.folder_dict["NumericalSim"], "3bp versus full perturbation model")

PF.Simple3DPlot(error_3bp[:,0], error_3bp[:,1], error_3bp[:,2], 
                title='Error vector over entire propagation',
                savefig = False, figfolder = '',
                name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]',
                zlabel='z [km]', xmod = 1000, ymod = 1000, zmod = 1000,
                xlim = None, ylim = None , zlim = None, fig = None, ax = None)

#%%

error_pos = np.linalg.norm(error_3bp,axis = 1)
error_pos2 = np.linalg.norm(error_fullmodel*EM,axis = 1)
error_baseline = np.linalg.norm(error_baseline_fullmodel,axis = 1)

t_days = np.linspace(t0,t1/(24*3600),error_3bp.shape[0])
fig,ax = PF.MakeFig(figsize=(12,3))
pp.plot(t_days[:250], error_pos2[:250], label = 'position error, S_1')
pp.plot(t_days[:250], error_pos[:250], label = 'position error, S_2',ls='--')
pp.plot(t_days[:250], error_baseline[:250], label = 'baseline error')
PF.DressFig(fig,ax, title = 'Satellite S1 error in barycentric frame', legendlist = ["x", "y", "z"],
            xlabel='Time in orbit [days]', ylabel = 'Error [m]',logScale = True)
pp.legend()
pp.tight_layout()
PF.saveFig(LL.folder_dict["NumericalSim"], "3pb versus full barycentric pos errors")

# fig,ax = PF.MakeFig(figsize=(12,3))
# PF.Matrix2DPlot(t_days[:50], np.abs(error_3bp[:,:50]), xmod = 1, ymod = 1000, fig=fig, ax = ax)
# PF.DressFig(fig,ax, title = 'Satellite S1 error in barycentric frame', legendlist = ["x", "y", "z"],
#             xlabel='Time in orbit [days]', ylabel = 'Error [km]',logScale = True)
# pp.tight_layout()
# PF.saveFig(LL.folder_dict["NumericalSim"], "3pb versus full barycentric pos errors focused")