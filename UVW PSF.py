# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 20:23:27 2020

@author: mrfan
"""

import numpy as np;
import matplotlib.pyplot as pp
from mpl_toolkits.mplot3d import Axes3D
import cv2

import os


def NormalizeArrayToDtype(arr, dType):
    
    im = np.array (arr/np.max(arr) * 255 , dtype = dType)
    return im

def ImagingAnalysis(uv_dist, I_sky):
    I_sky += + 1 *(np.abs(I_sky) == 0).astype(np.uint8) 
    fig = pp.figure(figsize=(4,4))
    pp.title('Original sky image')
    pp.axis([-150,150,-150,150])
    pp.imshow(I_sky,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    
    pp.savefig(figfolder + 'I_sky'+ figadd + '.png')
    
    L = np.floor(I_sky.shape[0]/2)
    I_fft = np.fft.fft2(I_sky)
    
    fft_shift = np.fft.fftshift(I_fft)
    im_fft = np.abs(fft_shift)
    im_fft = np.array( (im_fft)/np.max(im_fft) * 255 , dtype=np.uint8)

    fig = pp.figure(figsize=(4,4))
    pp.title('Centralized fourier transform of sky image')
    pp.axis([-150,150,-150,150])
    pp.imshow(NormalizeArrayToDtype(np.abs(im_fft),np.uint8),origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    pp.savefig(figfolder + 'CentralFourier'+ figadd +'.png')
    
    V = np.zeros(I_sky.shape,dtype=np.complex128)
    for i in range(0,uv_dist.shape[1]):
        u = np.int(uv_dist[0,i]/(1) + L) 
        v = np.int(uv_dist[1,i]/(1) + L)
        
        # V[u,v] = I_fft[u,v]
        V[v,u] = fft_shift[u,v]
        
        
    
    PSF = np.zeros(V.shape) + 1.* (np.abs(V) > 0) 
    I_PSF = NormalizeArrayToDtype(np.abs(
        np.fft.ifftshift(np.fft.ifft2(
            (PSF.astype(np.complex128))
                      ))),np.uint8)
    
    fig = pp.figure(figsize=(4,4))
    pp.title('Point Spread Function')
    pp.axis([-150,150,-150,150])
    pp.imshow(I_PSF,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    pp.savefig(figfolder + 'PSF'+ figadd +'.png')
    
    I_im = np.abs(np.fft.ifft2(np.fft.ifftshift(V)))
    
    fig = pp.figure(figsize=(4,4))
    pp.title('Sample function')
    pp.axis([-150,150,-150,150])
    pp.imshow(PSF,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    pp.savefig(figfolder + 'SampleFunc'+ figadd +'.png')
    
    fig = pp.figure(figsize=(4,4))
    pp.title('"Dirty" sky image')
    pp.axis([-150,150,-150,150])
    pp.imshow(NormalizeArrayToDtype(I_im, np.uint8),origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    pp.savefig(figfolder + 'I_obs'+ figadd +'.png')
    
    fig = pp.figure(figsize=(4,4))
    pp.title('Centralized fourier of dirty image')
    pp.axis([-150,150,-150,150])
    pp.imshow(NormalizeArrayToDtype(np.fft.fftshift(np.fft.fft2(I_im)), np.uint8),origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    pp.savefig(figfolder + 'I_fourier'+ figadd +'.png')
    
    
    C = NormalizeArrayToDtype(np.abs(I_im - I_PSF), np.uint8) 
    C[C<0] = 0
    fig = pp.figure(figsize=(4,4))
    pp.title('"Clean" sky image')
    pp.axis([-150,150,-150,150])
    pp.imshow(NormalizeArrayToDtype(C, np.uint8),origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    pp.savefig(figfolder + 'I_clean'+ figadd +'.png')
    
    
    
    return im_fft, V, I_im, I_PSF, C
    
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
    
    Projmatrix = np.zeros([2,data.shape[1]])
    
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
figfolder = cwd + '/Figures/uvw PSF/'
figadd = '_experiment'


if not os.path.exists(figfolder):
    os.makedirs(figfolder)
    print('made folder')
    
N= 30
    
sats = 150*(np.random.random([3,N])-0.5) 

# Spherical dist
# theta = 2*np.pi*np.random.random([2,N])
# sats = np.zeros((3,N))
# sats[0,:] = np.cos(theta[0])*np.sin(theta[1])
# sats[1,:] = np.sin(theta[0])*np.sin(theta[1])
# sats[2,:] = np.cos(theta[1])
# sats *= 50

bslns = GetBaselines(sats)
uvw_t = bslns

# Principle direction vectors
nx = np.array([0.0,1.0,0.0])
ny = np.array([-1.0,0.0,0.0])
nz = np.array([0.0,0.0,-1.0])

nw = np.array([1,1,1])

I_ref = cv2.imread(figfolder + 'SampleSky.jpg',cv2.IMREAD_GRAYSCALE)
projx = Project2Dplane(ny, uvw_t)

im_fft, V , I_im, I_psf , C= ImagingAnalysis(projx, I_ref)

c = np.int(uvw_t.shape[1]/2)



#%%

fig = pp.figure(figsize=(4,4))
pp.scatter(projx[0,0:c] , projx[1,0:c], c='b')
pp.scatter(projx[0,c:] , projx[1,c:], c='r')

if c < 50:
    for i in range(0,c):
#    pp.plot([projx[0,i], projx[0,c+i]], [projx[1,i],projx[1,c+i]], c='g', linestyle=':')
        pp.plot([0, projx[0,i]], [0,projx[1,i]], c='b', linestyle=':')
        pp.plot([0, projx[0,c+i]], [0, projx[1,c+i]], c='r', linestyle=':')


pp.title('uv Plane coverage')
pp.xlabel('x  [km]')
pp.ylabel('z  [km]')
pp.axis('equal')
pp.savefig(figfolder + 'uvwxproj'+ figadd + '.png')


fig = pp.figure(figsize=(4,4))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(sats[0,:], sats[1,:], sats[2,:])
ax.set_zlim([-50,50])
