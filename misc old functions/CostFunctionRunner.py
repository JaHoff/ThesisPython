# -*- coding: utf-8 -*-
"""
DEPRECATE - PART OF SCRAPPED PART OF THESIS

Minimal cost function runner file
Created on Tue Jun  2 17:15:28 2020

@author: Jurriez
"""

from scipy.special import j1 as bessel
from os import getcwd
import random

def NormalizeArrayToDtype(arr, dType):
    
    warn.filterwarnings("ignore")
    im = np.array(arr/np.max(arr) * 255 , dtype = dType)
    return im

def costFunctionAiry(psf):
    """Simple cost function for coarse PSF grids
    ripped from GF on 02-06"""
    W,H = psf.shape
    
    ang = np.degrees(np.arcsin(30/100e3))/2 #deg / px
    
    angles = np.empty_like(psf)
    x,y = np.arange(-W/2*ang,W/2*ang,ang),np.arange(-H/2*ang,H/2*ang,ang)
    x = x.reshape(1,len(x))
    y = y.reshape(len(y),1)
    
    WL = 30
    R = 100e3
    
    angles = np.sqrt(x**2 + y**2)
    x = 2*np.pi/WL * R * np.sin(np.deg2rad(angles))
    
    
    airy = (255*( 2* bessel(x)/x)**2).astype(int)
    cost = np.sum(np.sqrt((airy - psf)**2))/(W*H)
    
    return cost

def PSFanalysis(normals, samplepoints, resolution = 100., size=100e3):
    """Do a analysis of PSF for a series of view directions represented by a array of normal vectors \n
    Returns a complex128 poins spread function, and the uv sampling grid
    ripped from GF on 02-06"""
    # 
    
    L = int(size/resolution)    
    N = normals.shape[1]
    cost = np.zeros(N)
    
    grid, I_psf = np.zeros((2*L,2*L,N)),np.zeros((2*L,2*L,N))
    
    for i in range(0,N):
        # grid = np.zeros((2*L,2*L))
        Projection = Project2Dplane(normals[:,i], samplepoints)
        
        
        u = Projection[0,:]/(resolution)+L  # Was offset by +L to centralize the uv plane axis
        v = Projection[1,:]/(resolution)+L
            
        # Discard invalid baselines:
        v[ u < 0] = 5*L
        u[ u < 0] = 5*L
        u[ v < 0] = 5*L
        v[ v < 0] = 5*L
        ua, va = np.abs(u) ,np.abs(v)
        u[va >= 2*L] = 5*L
        v[va >= 2*L] = 5*L
        v[ua >= 2*L] = 5*L
        u[ua >= 2*L] = 5*L
            
        grid[v[v!=5*L].astype(int),u[u!=5*L].astype(int),i] = 1
 
        # Clear up as much memory as possible for the next step
        del v,u, Projection
        
        I_psf[:,:,i] = np.fft.ifft2(grid[:,:,i])
        
        cost[i] = costFunctionAiry(NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_psf[:,:,i])), np.uint8))
        print("Now at normal index " + str(i) + " of " + str(N) + " yielded a cost of " + str(cost[i]))
    
    return I_psf,grid,cost

def Project2Dplane(normal, data):
    """Project a series of three dimensional points onto a plane orthogonal to the normal vector \n
        IN \n
        normal: a normalized 3 dimensional direction vector  towards the viewing plane \n
        data: Collection of 3-dimensional sample points to be projected \n
        OUT \n
        Projmatrix: 2 dimensional projected matrix of baseline data. \n
    """

#    Normalize the direction vectors, also define the mystery direction that is the second axis.
    normal = normal/np.linalg.norm(normal)
    if (np.abs(normal[2]) == 1.0): 
        normal_z = np.array([1.0,0.0,0.0])
    else:
        normal_z = np.array([0.0,0.0,1.0]) 
        
    # Turn the z-direction vector into an orthogonal vector to the normal

    normal_z -= normal_z.dot(normal)*normal
    normal_z /= np.linalg.norm(normal_z)
    # Finish the base axes by adding the final vector
    normal_e = np.cross(normal,normal_z)
    
    Projmatrix = np.zeros([2,data.shape[0]])
    
    # There might be a way to do this in a single matrix computation
    for i in range(0,data.shape[0]):
        q_proj = data[i,:] - np.dot(data[i,:], normal) * normal # make vectors orthogonal to n
        Projmatrix[1,i] = np.dot(normal_e,q_proj)
        Projmatrix[0,i] = np.dot(normal_z,q_proj)    
    return Projmatrix



def sampleSphereFib(N, randomize = False):
    """Create a set of points along a sphere, distributed following a Fibbonaci lattice \n
    Returns a set of 3-dimensional cartesian points."""
    
    rnd = 1.
    if randomize:
        rnd = random.random() * N

    points = []
    offset = 2./N
    increment = np.pi * (3. - np.sqrt(5.));
    
    points = np.zeros((3,N))

    for i in range(N):
        y = ((i * offset) - 1) + (offset / 2);
        r = np.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % N) * increment

        x = np.cos(phi) * r
        z = np.sin(phi) * r

        points[0,i] = x
        points[1,i] = y
        points[2,i] = z
        
    return points


import numpy as np 

# cwd = getcwd()
# datafolder = cwd + '/Data/PerturbationAnalysis/'
# File = datafolder + 'propagationHistory_0_fullprop.dat'

# print("Importing data")
# t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,1, True)
# print("Imported data succesfully")

# #%%
# sat = np.concatenate((x,y,z), axis=1)
# co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)

# N = co_x.shape[0]

# #%% 
# BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)

# normals = GF.sampleSphereFib(30)

# I_PSF, I_S,cost = GF.PSFanalysis(normals, BL_c,resolution = res)

# cost_total = np.sum(cost)