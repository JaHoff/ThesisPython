# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:27:28 2019

@author: USER
"""

import numpy as np;
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os



def Project2Dplane(normal, data):
    
#    Normalize the direction vectors
    normal = normal/np.linalg.norm(normal)
    if (np.abs(normal[2]) == 1.0): 
        normal_z = np.array([1.0,0.0,0.0])
    else:
        normal_z = np.array([0.0,0.0,1.0])
        
    # Turn the z-direction vector into an orthogonal vector to the normal
    normal_z -= normal_z.dot(normal)*normal
    normal_z /= np.linalg.norm(normal_z)
    # Finish the base axes by adding the final vector
    normal_e =np.cross(normal,normal_z)
    
    Projmatrix = np.zeros([3,data.shape[1]])
    
    # There might be a way to do this in a single matrix computation
    for i in range(0,data.shape[1]):
        q_proj = data[:,i] - np.dot(data[:,i], normal) * normal # make vectors orthogonal to n
        Projmatrix[0,i] = np.dot(normal_e,q_proj)
        Projmatrix[1,i] = np.dot(normal_z,q_proj)    
    return Projmatrix

def GetBaselines(sats):
    # sats is a 3xN array of 3-dimensional satellite 
    N = sats.shape[1] #5
    
    baselines = np.zeros((3,N*(N-1)))
    n = np.int32(N*(N-1)/2)
    for i in range(0,N-1):
        s = i*N - np.int32( (i**2 + i)/2 )
        e = (i+1)*N - np.int32(((i+1)**2 +(i+1))/2)
        baselines[:,s : e] = sats[:,(i+1):] - (np.ones((N-(i+1),1))*sats[:,i]).T
        
    baselines[:, n:] = -baselines[:,:n]
    return baselines

cwd = os.getcwd()
figfolder = cwd + '/Figures/uvw_geometry/'


if not os.path.exists(figfolder):
    os.makedirs(figfolder)
    print('made folder')


sats = np.random.normal(0,100,[3,5])

bslns = GetBaselines(sats)
uvw_t = bslns


nx = np.array([0.0,1.0,0.0])
ny = np.array([-1.0,0.0,0.0])
nz = np.array([0.0,0.0,-1.0])

nw = np.array([1,1,1])


projx = Project2Dplane(nx, uvw_t)
projy = Project2Dplane(ny, uvw_t)
projz = Project2Dplane(nz, uvw_t)
projw = Project2Dplane(nw,uvw_t)


c = np.int(uvw_t.shape[1]/2)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(uvw_t[0,:c], uvw_t[1,:c], uvw_t[2,:c], s=20, c='b', depthshade=True, label = 'uvw baselines')
ax.scatter(uvw_t[0,c:], uvw_t[1,c:], uvw_t[2,c:], s=20, c='r', depthshade=True, label = 'negative pairs')

for i in range(0,c):
    ax.plot([uvw_t[0,i], uvw_t[0,c+i]], [uvw_t[1,i],uvw_t[1,c+i]], [uvw_t[2,i], uvw_t[2,c+i]], c='g', linestyle=':')
    
ax.set_title('uvw space baselines')    
ax.set_xlabel('x [km]')
ax.set_ylabel('y [km]')
ax.set_zlabel('z [km]')

plt.savefig(figfolder + 'uvw3d.png')


fig = plt.figure(2)
plt.scatter(projx[0,0:c] , projx[1,0:c], c='b')
plt.scatter(projx[0,c:] , projx[1,c:], c='r')

for i in range(0,c):
#    plt.plot([projx[0,i], projx[0,c+i]], [projx[1,i],projx[1,c+i]], c='g', linestyle=':')
    plt.plot([0, projx[0,i]], [0,projx[1,i]], c='b', linestyle=':')
    plt.plot([0, projx[0,c+i]], [0, projx[1,c+i]], c='r', linestyle=':')

plt.title('xz-plane projection')
plt.xlabel('x  [km]')
plt.ylabel('z  [km]')

plt.savefig(figfolder + 'uvwxproj.png')

fig = plt.figure(3)
plt.scatter(projy[0,0:c] , projy[1,0:c], c='b')
plt.scatter(projy[0,c:] , projy[1,c:], c='r')

for i in range(0,c):
#    plt.plot([projy[0,i], projy[0,c+i]], [projy[1,i],projy[1,c+i]], c='g', linestyle=':')
    plt.plot([0, projy[0,i]], [0,projy[1,i]], c='b', linestyle=':')
    plt.plot([0, projy[0,c+i]], [0, projy[1,c+i]], c='r', linestyle=':')
    
plt.title('yz-plane projection')
plt.xlabel('y  [km]')
plt.ylabel('z  [km]')
plt.savefig(figfolder + 'uvwyproj.png')

fig = plt.figure(4)
plt.scatter(projz[0,0:c] , projz[1,0:c], c='b')
plt.scatter(projz[0,c:] , projz[1,c:], c='r')

for i in range(0,c):
#    plt.plot([projz[0,i], projz[0,c+i]], [projz[1,i],projz[1,c+i]], c='g', linestyle=':')
    plt.plot([0, projz[0,i]], [0,projz[1,i]], c='b', linestyle=':')
    plt.plot([0, projz[0,c+i]], [0, projz[1,c+i]], c='r', linestyle=':')
    
plt.title('yx-plane projection')
plt.xlabel('y  [km]')
plt.ylabel('x  [km]')
plt.savefig(figfolder + 'uvwzproj.png')


fig = plt.figure(5, figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(sats[0,:] , sats[1,:], sats[2,:], c='k', label='Satellites')

for i in range(0,4):
    for j in range(4,i,-1):
        line,= ax.plot([sats[0,i], sats[0,j]], [sats[1,i],sats[1,j]], [sats[2,i], sats[2,j]], c='g', linestyle=':')

line.set_label('Physical baseline')
plt.title('Satellite constellation')
plt.xlabel('y [km]')
plt.ylabel('x [km]')
ax.set_zlabel('z [km]')
plt.legend()
plt.savefig(figfolder + 'uvwsatgeo.png')