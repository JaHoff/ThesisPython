# -*- coding: utf-8 -*-
"""
Experiment with the gradient of the potential field theta
Created on Tue Apr 28 15:53:16 2020

@author: Jurriez
"""

import numpy as np
import sys
from os import chdir
from os import getcwd
from matplotlib import cm
import matplotlib.pyplot as pp

sys.path.append('../')
chdir('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF
import MayaviPlotting as MF
chdir('Potential Gradient')

#%% Theta and gradient functions
def theta(x,y,z):
    
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = 1- m_E/(m_t)

    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    
    r1 = np.sqrt( (x+mu)**2 + y**2 + z**2)
    r2 = np.sqrt( (x-1+mu)**2 + y**2 + z**2)
    
    thetamap = (1-mu)/r1 + mu/r2 + x**2/2 + y**2/2
    
    return thetamap

def dthetadx(x,y,z):
    
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = 1- m_E/(m_t)

    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    
    thetamap = x -( 1-mu)*(2*x+2*mu)/(2*((x+mu)**2+y**2+z**2)**(3/2))-mu*(2*x-2+2*mu)/(2*((x-1+mu)**2+y**2+z**2)**(3/2))
    
    return thetamap

def dthetady(x,y,z):
    
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = 1- m_E/(m_t)

    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    
    thetamap = y -(1-mu)*y/((x+mu)**2+y**2+z**2)**(3/2)-mu*y/((x-1+mu)**2+y**2+z**2)**(3/2)
    
    return thetamap

def dthetadz(x,y,z):
    
    
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = 1- m_E/(m_t)

    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    
    thetamap = -(1-mu)*z/((x+mu)**2+y**2+z**2)**(3/2)-mu*z/((x-1+mu)**2+y**2+z**2)**(3/2)
    
    return thetamap

def dthetadxdy(x,y,z):
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = 1- m_E/(m_t)
    
    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    
    thetamap = 3*(1 - mu)*(2*x + 2*mu)*y/(2*((x + mu)**2 + y**2 + z**2)**(5/2)) + (3*mu*(2*x - 2 + 2*mu)*y)/(2*((x - 1 + mu)**2 + y**2 + z**2)**(5/2))

    return thetamap

#%% compute functions
def dthetadp(x,y,z):
    
    dtdx = dthetadx(x,y,z)
    dtdy = dthetady(x,y,z)
    dtdz = dthetadz(x,y,z)
    
    dtheta = np.sqrt( dtdx**2 + dtdy**2 + dtdz**2)    
    
    
    #dtheta = dtdx + dtdy + dtdz
    return dtheta, dtdx, dtdy, dtdz

def computeGradientToPoint(x,y,z,xt,yt,zt,dtdx,dtdy,dtdz):
    # Gradient approach
    xd,yd,zd = xt*np.ones(x.shape), yt*np.ones(y.shape), zt*np.ones(z.shape)
    
    x = xd - x
    y = yd - y
    z = zd - z
    
    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    r = np.sqrt(x**2 + y**2 + z**2)
    dtheta = dtdx*x/r + dtdy*y/r + dtdz*z/r
    
    return dtheta

def computeVCPF(dtdx, dtdy,dtdz):
    G = 6.674e-11
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = m_E/(m_t)
    
    EM = 364522821.08369416
    
    tc = np.sqrt(EM**3/(G*m_t))
    
    
    dtdxp = 2*dtdy + dtdx

    dtdyp = -2*dtdx+ dtdy
    
    dtdzp = dtdz

    vcpff = np.sqrt(dtdxp**2 + dtdyp**2 + dtdzp**2)
    
    return dtdxp, dtdyp, dtdzp, vcpff

def computeVCPFSz(dtdx, dtdy,dtdz,Sz):
    G = 6.674e-11
    m_E = 5.97237e24
    m_M = 7.34767309e22
    m_t = m_E + m_M
    mu = m_E/(m_t)
    
    EM = 364522821.08369416
    
    tc = np.sqrt(EM**3/(G*m_t))
    
    scale = np.sqrt(2*Sz/np.sqrt(dtdx**2 + dtdy**2 + dtdz**2)) 
    
    dtdxp = 2*dtdy*scale+ dtdx

    dtdyp = -2*dtdx*scale+ dtdy
    
    dtdzp = dtdz

    vcpff = np.sqrt(dtdxp**2 + dtdyp**2 + dtdzp**2)
    
    return dtdxp, dtdyp, dtdzp, vcpff

def computeRelVectorProj(x0,y0,z0, x,y,z, u,v,w):
    x = x.reshape(len(x),1,1)
    y = y.reshape(1,len(y),1)
    z = z.reshape(1,1,len(z))
    
    xr = x-x0
    yr = y - y0
    zr = z - z0
    
    rb = np.sqrt(xr**2 + yr**2 + zr**2)
    
    dotproduct = (xr*u + yr*v + zr*w)/rb
    
    return dotproduct
    
    
MF.CloseAll()
#%% Set up solution space
# test_mesh()

figuresize = (1800,1200)
x,y,z = (np.arange(-1.2,1.21,0.01) for i in range(0,3))
N = len(z)
N2 = int(N/2)

#%% set up mu

m_E = 5.97237e24
m_M = 7.34767309e22
m_t = m_E + m_M
mu = 1- m_E/(m_t)

#%%

thetamap = theta(x,y,z)

dthetamap, dtdx, dtdy, dtdz = dthetadp(x,y,z)




# #%%
# dtdxdy = dthetadxdy(x,y,z)
# xp = np.outer(x,np.ones(len(y)))
# yp = np.outer(y,np.ones(len(x))).T
# zp = np.abs(dtdxdy[:,:,N2])
# zp[zp>3] = 3
# zp[zp<-3] = -3

# fig = MF.MakeFig("Omega dx dy diff",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
# MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
# zp += 0.01
# MF.PlotColorSphere(fig, 1, x=0,y=0,z=1,R = 0.1, c = (0,1,0), opacity = 1)
# MF.PlotColorSphere(fig, 1, x=1,y=0,z=1,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# # plot lagrangian points
# l4_scale = 0.01
# l4_color = (1,0,0)
# MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=1,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=1,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=1,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=1,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=1,R = l4_scale, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
# MF.Title(fig,"Omega differential in x and y")


#%% Mayavi surface plot - potential gradient quivers
xp = np.outer(x,np.ones(len(y)))
yp = np.outer(y,np.ones(len(x))).T
zp = 2-thetamap[:,:,N2]
zp[zp < -0] = -0 # 0.5 / 2.5 makes for a nice solution space manifold

elementZ = np.max(zp)

# Establish a velocity corrected potential field

dtdxp, dtdyp, dtdzp, vcpff = computeVCPF(dtdx,dtdy, dtdz)
relL4Velocity = computeRelVectorProj(0.5-mu,np.sqrt(3)/2,0, x,y,z, dtdxp,dtdyp,dtdzp) # used for rel velocity to L4

# curb the magnitude for the figure
lim = 1.5
dtdxp[dtdxp < -lim] = -lim
dtdxp[dtdxp > lim] = lim
dtdyp[dtdyp < -lim] = -lim
dtdyp[dtdyp > lim] = lim
dtdzp[dtdzp < -lim] = -lim
dtdzp[dtdzp > lim] = lim

fig = MF.MakeFig("Potential vector field",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
zp += 0.01
# MF.Plot3DQuiver(fig,dtdxp[:,:,N2],dtdyp[:,:,N2],dtdzp[:,:,N2],xp,yp,zp,sf = 0.05)
MF.PlotColorSphere(fig, 1, x=0-mu,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1-mu,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404-mu,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596-mu,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071-mu,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field gradient with VCPF vectors")

#%% Mayavi surface plot - stable velocity quivers
zp = 2-thetamap[:,:,N2]
zp[zp < -0] = -0 # 0.5 / 2.5 makes for a nice solution space manifold

elementZ = np.max(zp)

dtdxp = np.copy(dtdx)
dtdyp = np.copy(dtdy)
dtdzp = np.copy(dtdz)

tester = computeRelVectorProj(0.5-mu,np.sqrt(3)/2,0, x,y,z, 0.5*dtdy,-0.5*dtdx,dtdz)
tester[tester > 0.4] = 0.4
tester[tester < -.4] = -.4
fig = MF.MakeFig("Rel. stable velocity from L4", bgcolor=(1,1,1), fgcolor=(0,0,0))
MF.PlotSurface(fig,xp,yp,tester[:,:,N2], cmap = "copper")
MF.AddDefaultAxes(fig)


# curb the magnitude for the figure
lim = 1.5
dtdxp[dtdxp < -lim] = -lim
dtdxp[dtdxp > lim] = lim
dtdyp[dtdyp < -lim] = -lim
dtdyp[dtdyp > lim] = lim
dtdzp[dtdzp < -lim] = -lim
dtdzp[dtdzp > lim] = lim

fig = MF.MakeFig("Potential field - stable vel",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
zp += 0.01
MF.Plot3DQuiver(fig,0.5*dtdyp[:,:,N2],-0.5*dtdxp[:,:,N2],dtdzp[:,:,N2],xp,yp,zp,sf = 1)
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field gradient with stable velocity vectors")

Stable_velocity_magnitude = np.sqrt((0.5*dtdyp)**2 + (-0.5*dtdxp)**2 + dtdzp**2)


#%% Mayavi surface plot - VCPF field
zp = 2-vcpff[:,:,N2]
zp[zp < -0] = -0 # 0.5 / 2.5 makes for a nice solution space manifold

fig = MF.MakeFig("Potential field - VCPFF",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"VCPF Field")



#%% Mayavi surface plot - stable velocity magnitude
# xp = np.outer(x,np.ones(len(y)))
# yp = np.outer(y,np.ones(len(x))).T
# zp = Stable_velocity_magnitude[:,:,N2]

# elementZ = np.max(zp)


# fig = MF.MakeFig("Stable velocity magnitude",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
# MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
# zp += 0.01
# MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
# MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# # plot lagrangian points
# l4_scale = 0.01
# l4_color = (1,0,0)
# MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
# MF.Title(fig,"Magnitude of the stable velocity field")
#%% Plot VCPF minus stable velocity vectors

fig = MF.MakeFig("Potential field - VCPF minus stable velocity vectors",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
zp += 0.01
MF.Plot3DQuiver(fig,dtdxp[:,:,N2]-0.5*dtdyp[:,:,N2],dtdyp[:,:,N2]+0.5*dtdxp[:,:,N2],np.zeros(dtdxp[:,:,N2].shape),xp,yp,zp,sf = 0.05)
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field gradient with VCPF vectors minus stable velocities")

#%% VCPFF Field figure
zp = 2-thetamap[:,:,N2]
cutoff = 0.05
zp[vcpff[:,:,N2] > cutoff] = 0 # 0.5 / 2.5 makes for a nice solution space manifold

fig = MF.MakeFig("VCPFF corrected area",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,f"Potential field where VCPFF < {cutoff}")



#%% Compare gradient towards L4 to that towards Earth and moon

grad_L4 = computeGradientToPoint(x,y,z,0.5-mu,np.sqrt(3)/2,0,dtdx,dtdy,dtdz)

r = np.sqrt(x**2 + y**2 + z**2)
grad_Earth = computeGradientToPoint(x,y,z,-mu,0,0,dtdx,dtdy,dtdz)

xm = 1-x
r = np.sqrt(xm**2 + y**2 + z**2)
grad_Moon = computeGradientToPoint(x,y,z,1-mu,0,0,dtdx,dtdy,dtdz)

grad_L4 = grad_L4[:,:,N2]
grad_Earth = grad_Earth[:,:,N2]
grad_Moon = grad_Moon[:,:,N2]



#%% Mayavi surface plot
# Alter map to only show spaces with potential field smaller than 0.5
zp = -grad_Earth
zp[zp < -0.5] = -0.5 # 0.5 / 2.5 makes for a nice solution space manifold
zp[zp>0.5] = 0.5
elementZ = 0

fig = MF.MakeFig("Potential field - Earth",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field gradient towards Earth - inverted")

#%% Mayavi surface plot
# Alter map to only show spaces with potential field smaller than 0.5
zp = -grad_Moon
zp[zp < -0.5] = -0.5 # 0.5 / 2.5 makes for a nice solution space manifold
zp[zp>0.5] = 0.5
elementZ = 0

fig = MF.MakeFig("Potential field - Moon",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field gradient towards Moon - inverted")

#%% Mayavi surface plot
# Alter map to only show spaces with potential field smaller than 0.5
zp = -grad_L4
zp[zp < -0.5] = -0.5 # 0.5 / 2.5 makes for a nice solution space manifold
zp[zp>0] = 0
elementZ = 0

fig = MF.MakeFig("Potential field - L4",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field gradient towards L4")

#%% Alter L4 to exclude false positives

grad_L4[np.abs(grad_L4) < np.abs(grad_Earth)] = -10
grad_L4[np.abs(grad_L4) < np.abs(grad_Moon)] = -10
grad_L4[yp < 0] = -10


#%% Mayavi surface plot
# Alter map to only show spaces with potential field smaller than 0.5
zp = 2-thetamap[:,:,N2]
zp[zp < 0] = 0 # 0.5 / 2.5 makes for a nice solution space manifold
elementZ = np.max(zp)

fig = MF.MakeFig("Potential field - inverted",bgcolor = (1,1,1) , fgcolor = (0,0,0),size=figuresize)
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field - inverted: 2 - \Omega")

#%% Mayavi surface plot
zp = thetamap[:,:,N2]-1
zp[zp > 0.5] = 0
elementZ = np.max(zp)

fig = MF.MakeFig("Potential field - components smaller than 1.5, z = 0",bgcolor = (1,1,1) , fgcolor = (0,0,0))
MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig,"Potential field - components smaller than 1.5, z = 0")

#%%
dtdxa = np.abs(dtdx)
dtdya = np.abs(dtdy)
dtdza = np.abs(dtdz)
test = np.zeros(dtdx.shape)
test[dtdxa > test] = dtdxa[dtdxa > test]
test[dtdya > test] = dtdya[dtdya > test]
test[dtdza > test] = dtdza[dtdza > test]

#%% Mayavi surface plot -  testL4 gradient dominance

# zp = 1-test[:,:,N2]
# zp[zp < 0.95] = 0.5
# zp[zp > 0.95] = 1

# fig = MF.MakeFig("Gradient approx 0" ,bgcolor = (1,1,1) , fgcolor = (0,0,0))
# MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')
# MF.PlotColorSphere(fig, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
# MF.PlotColorSphere(fig, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# # plot lagrangian points
# l4_scale = 0.01
# l4_color = (1,0,0)
# MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=0.5,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig, 1, x=0.5,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
# MF.Title(fig,"points where maximum gradient component \approx 0")

#%% Mayavi surface plot - potential gradient
zp = grad_L4
zp[zp < -0.10] = -.10
zp[zp > 0.1] = 0.1



elementZ = np.max(zp)

fig2 = MF.MakeFig('Potential gradient pull towards L4, capped at 0.1',bgcolor = (1,1,1) , fgcolor = (0,0,0))
MF.PlotSurface(fig2, xp,yp,zp, lw = 4.0, cmap = 'copper')
MF.PlotColorSphere(fig2, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig2, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig2, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig2, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig2, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
MF.Title(fig2,'Potential gradient pull towards L4, capped at 0.1')

# #%% Mayavi surface plot - potential gradient dx
# zp = dtdx[:,:,N2]
# zp[zp <-1] = -1
# zp[zp > 1] = 1


# elementZ = np.max(zp)

# fig2 = MF.MakeFig('dtdx gradient',bgcolor = (1,1,1) , fgcolor = (0,0,0))
# MF.PlotSurface(fig2, xp,yp,zp, lw = 4.0, cmap = 'copper')
# MF.PlotColorSphere(fig2, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# # plot lagrangian points
# l4_scale = 0.01
# l4_color = (1,0,0)
# MF.PlotColorSphere(fig2, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)

#%% Mayavi surface plot - potential gradient dy
# zp = dtdx[:,:,N2] + dtdy[:,:,N2]
# zp[zp <-1] = -1
# zp[zp > 1] = 1

# elementZ = np.max(zp)

# fig2 = MF.MakeFig('dtdy + dtdx gradient',bgcolor = (1,1,1) , fgcolor = (0,0,0))
# MF.PlotSurface(fig2, xp,yp,zp, lw = 4.0, cmap = 'copper')
# MF.PlotColorSphere(fig2, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# # plot lagrangian points
# l4_scale = 0.01
# l4_color = (1,0,0)
# MF.PlotColorSphere(fig2, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)
# MF.Title(fig2,"X and Y gradients addition")

#%% Mayavi surface plot - potential gradient dy
# zp = dtdy[:,:,N2]
# zp[zp <-1] = -1
# zp[zp > 1] = 1

# elementZ = np.max(zp)

# fig2 = MF.MakeFig('dtdy gradient',bgcolor = (1,1,1) , fgcolor = (0,0,0))
# MF.PlotSurface(fig2, xp,yp,zp, lw = 4.0, cmap = 'copper')
# MF.PlotColorSphere(fig2, 1, x=0,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=1,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# # plot lagrangian points
# l4_scale = 0.01
# l4_color = (1,0,0)
# MF.PlotColorSphere(fig2, 1, x=0.8404,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=1.1596,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=-1.0071,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
# MF.PlotColorSphere(fig2, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)



#%% Close up of the L4 point


d = 0.05
dt = d/500
x = np.arange(0.5-mu-d,0.5-mu+d,dt)
y = np.arange(np.sqrt(3)/2-d,np.sqrt(3)/2+d, dt )
z = np.array([-d,0,d])

Nx = len(x)
N = len(z)
N2 = int(N/2)

dthetamapcu,dx,dy,dz = dthetadp(x,y,z)





zp = computeGradientToPoint(x,y,z,0.5-mu,np.sqrt(3)/2,0,dx,dy,dz)[:,:,1]

x -= 0.5-mu
y -= np.sqrt(3)/2

EM = 364522821.08369416
xp = np.outer(x,np.ones(len(y))) * EM/1000
yp = np.outer(y,np.ones(len(x))).T * EM/1000

tester = computeRelVectorProj(0.5-mu,np.sqrt(3)/2,0, x,y,z, 0.5*dy,-0.5*dx,dz)
tester[np.abs(tester) > 2] = 2
fig = MF.MakeFig("Rel. stable velocity from L4 - close up", bgcolor=(1,1,1), fgcolor=(0,0,0))
MF.PlotSurface(fig,xp,yp,np.abs(tester[:,:,N2]), cmap = "copper",warp_scale = 'auto')
MF.AddDefaultAxes(fig)


#%%

# zp = -dthetamapcu[:,:,N2]
# zp[zp <-1] = -1
# zp[zp > 0.1] = 0.1

# zp *= 1e5
#%% Mayavi surface plot
# elementZ = np.max(zp)

# fig3 = MF.MakeFig('Close up of L4 gradient to L4',bgcolor = (1,1,1) , fgcolor = (0,0,0))
# MF.PlotSurface(fig3, xp,yp,zp, lw = 1.0, cmap = 'copper', warp_scale = 1)
# # Add L4
# MF.PlotColorSphere(fig3, 1, x=0,y=0,z=elementZ,R = 100, c = l4_color, opacity = 1)

# # Add axes

# MF.AddDefaultAxes(fig3, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)


#%% 

meshmap = 0*np.ones(thetamap.shape)
# meshmap[thetamap <= 1.5] = 1
meshmap[Stable_velocity_magnitude < 0.015] = 1
x,y,z = np.mgrid[-1.2:1.2:241j, -1.2:1.2:241j, -1.2:1.2:241j]

#%% 
fig = MF.MakeFig("Suitable solution manifold - stable velocity field", bgcolor=(1,1,1), fgcolor=(0,0,0))
MF.Plot3DContour(fig,meshmap, cmap = 'copper', lw = 2.0, opacity =0.8 ,x = x, y = y, z = z)
MF.PlotColorSphere(fig, 1, x=0,y=0,z=0,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=0,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.05
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)

#%%  relative velocity test



fig = MF.MakeFig("VCPF relative to L4",bgcolor=(1,1,1), fgcolor=(0,0,0),size=figuresize)
MF.PlotSurface(fig,x[:,:,N2], y[:,:,N2],np.abs(relL4Velocity[:,:,N2]))

MF.PlotColorSphere(fig, 1, x=0,y=0,z=0,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=0,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.05
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)

MF.AddDefaultAxes(fig)
#%% Make a selection based on the vcpff map
meshmap = 0*np.ones(thetamap.shape)
meshmap[vcpff < cutoff] = 1

fig = MF.MakeFig("Suitable solution manifold - VCPFF", bgcolor=(1,1,1), fgcolor=(0,0,0))
MF.Plot3DContour(fig,meshmap, cmap = 'copper', lw = 2.0, opacity =0.8 ,x = x, y = y, z = z)
MF.PlotColorSphere(fig, 1, x=0,y=0,z=0,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=0,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.05
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)

#%% Adjusted potential mnus distance

L4 = np.array([0.5-mu, np.sqrt(3)/2, 0])

Dx = x - L4[0]
Dy = y - L4[1]
Dz = z - L4[2]

adj_theta = thetamap + np.sqrt(Dx**2 + Dy**2 + Dz**2)





meshmap = 0*np.ones(thetamap.shape)
meshmap[adj_theta > 1.8] = 1

fig = MF.MakeFig("Suitable solution manifold - spatially adjusted theta", bgcolor=(1,1,1), fgcolor=(0,0,0))
MF.Plot3DContour(fig,meshmap, cmap = 'plasma', lw = 2.0, opacity =0.5 ,x = x, y = y, z = z)
MF.PlotColorSphere(fig, 1, x=0,y=0,z=0,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1,y=0,z=0,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.05
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071,y=0,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=0,R = l4_scale, c = l4_color, opacity = 1)

# Add axes

MF.AddDefaultAxes(fig2, color = (0,0,0), xlabel = 'x', ylabel = 'y', zlabel = 'z', Nlabels=5)

#%% Make semi-2d quiver plot
x,y,z = np.mgrid[-1.2:1.2:41j, -1.2:1.2:41j, -0:0:1j]


elementZ = np.max(z)

fig = MF.MakeFig("semi2DQuiver", figuresize, bgcolor = (1,1,1), fgcolor = (0,0,0))


dm, dtdx,dtdy,dtdz = dthetadp(x[:,0,0],y[0,:,0],z[0,0,:])
dtdxp = np.copy(dtdx)
dtdyp = np.copy(dtdy)
dtdzp = np.copy(dtdz)
# curb the magnitude for the figure
lim = 0.5
dtdxp[dtdxp < -lim] = -lim
dtdxp[dtdxp > lim] = lim
dtdyp[dtdyp < -lim] = -lim
dtdyp[dtdyp > lim] = lim
dtdzp[dtdzp < -lim] = -lim
dtdzp[dtdzp > lim] = lim
# MF.PlotSurface(fig, xp,yp,zp, lw = 4.0, cmap = 'copper')

MF.Plot3DQuiver(fig,0.5*dtdyp[:,:,0],-0.5*dtdxp[:,:,0],dtdzp[:,:,0],x[:,:,0],y[:,:,0],z[:,:,0],sf = 0.03,cmap = 'plasma')
MF.PlotColorSphere(fig, 1, x=0-mu,y=0,z=elementZ,R = 0.1, c = (0,1,0), opacity = 1)
MF.PlotColorSphere(fig, 1, x=1-mu,y=0,z=elementZ,R = 0.05, c = (0.8,0.8,0.8), opacity = 1)

# plot lagrangian points
l4_scale = 0.01
l4_color = (1,0,0)
MF.PlotColorSphere(fig, 1, x=0.8404-mu,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=1.1596-mu,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=-1.0071-mu,y=0,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)
MF.PlotColorSphere(fig, 1, x=0.5-mu,y=-np.sqrt(3)/2,z=elementZ,R = l4_scale, c = l4_color, opacity = 1)


MF.AddDefaultAxes(fig)
#%%
fig,ax = PF.MakeFig(figsize = (8,8))
ax.quiver(x[:,:,0], y[:,:,0], 0.5*dtdyp[:,:,0],-0.5*dtdxp[:,:,0], pivot = 'mid')
ax.scatter(0.5-mu, np.sqrt(3)/2, color='r')
ax.scatter(0.5-mu, -np.sqrt(3)/2, color='r')
ax.scatter(0.8404, 0, color='r')
ax.scatter(1.1596, 0, color='r')
ax.scatter(-1.00711, 0, color='r')
ax.scatter(1-mu, 0, color='b',s=120)
ax.scatter(-mu, 0, color='b', s=240)

PF.DressFig(fig,ax, 'Stable velocity field in Earth-Moon system', True, LL.folder_dict["Methodology"], name = 'Stable velocity field 3BP',
            xlabel = 'x', ylabel = 'y')

fig,ax = PF.MakeFig(figsize = (8,8))
ax.quiver(x[:,:,0], y[:,:,0], dtdxp[:,:,0],dtdyp[:,:,0], pivot = 'mid')
ax.scatter(0.5-mu, np.sqrt(3)/2, color='r')
ax.scatter(0.5-mu, -np.sqrt(3)/2, color='r')
ax.scatter(0.8404, 0, color='r')
ax.scatter(1.1596, 0, color='r')
ax.scatter(-1.00711, 0, color='r')
ax.scatter(1-mu, 0, color='b',s=120)
ax.scatter(-mu, 0, color='b', s=240)

PF.DressFig(fig,ax, 'Initial acceleration for motionless object', True, LL.folder_dict["Methodology"], name = 'Initial acceleration field 3BP',
            xlabel = 'x', ylabel = 'y')

dtdxp, dtdyp, dtdzp, vcpff = computeVCPFSz(dtdx,dtdy, dtdz, 0.06)

# curb the magnitude for the figure
lim = 0.5
dtdxp[dtdxp < -lim] = -lim
dtdxp[dtdxp > lim] = lim
dtdyp[dtdyp < -lim] = -lim
dtdyp[dtdyp > lim] = lim
dtdzp[dtdzp < -lim] = -lim
dtdzp[dtdzp > lim] = lim

fig,ax = PF.MakeFig(figsize = (8,8))
ax.quiver(x[:,:,0], y[:,:,0], dtdxp[:,:,0],dtdyp[:,:,0], pivot = 'mid')
ax.scatter(0.5-mu, np.sqrt(3)/2, color='r')
ax.scatter(0.5-mu, -np.sqrt(3)/2, color='r')
ax.scatter(0.8404, 0, color='r')
ax.scatter(1.1596, 0, color='r')
ax.scatter(-1.00711, 0, color='r')
ax.scatter(1-mu, 0, color='b',s=120)
ax.scatter(-mu, 0, color='b', s=240)

PF.DressFig(fig,ax, 'Auto-velocity acceleration in Earth-Moon system', True, LL.folder_dict["Methodology"], name = 'Natural acceleration 3BP',
            xlabel = 'x', ylabel = 'y')
# PF.Simple2DScatter(np.array([0.5-mu]), np.array([np.sqrt(3)/2]), color = 'r',fig = fig, ax = ax )

#%% Convert the existing gradient to a velocity gradient, an attempt based on what I think is right

# Problem: Scale of resultts doesn' match what one would expect...

dtdxp, dtdyp, dtdzp, vcpff = computeVCPFSz(dtdx,dtdy, dtdz, 0.06)

RE = LL.Constants["Lunar_semi_major_axis"]
G = LL.Constants["G"]
mD = 5.972E24 + 7.34767309E22

t_char = np.sqrt( (RE**3) /(G*mD))

accel_factor = RE / t_char**2


accelerations = dthetamap*accel_factor