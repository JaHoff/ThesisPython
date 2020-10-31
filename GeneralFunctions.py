# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:37:31 2020

@author: USER
"""
import numpy as np
import random
#import scipy.fftpack as fft

import warnings as warn

def GetBaseLines(x,y,z, vx, vy, vz, t):
    """Get the baselines from a series of points spread in cartesian three-dimensional space\n
    Returns baseline coordinates, velocities and magnitudes spread and in combined arrays: \n
    BLx, BLy, BLz, BLc, BLrx, BLry, BLrz, BLrc, blM, blrM"""
    Width = x.shape[1]
    Height = x.shape[0]
    n = np.int32(Width*(Width-1)/2)
    
    BLx, BLy, BLz, blM, BLrx, BLry, BLrz, blrM = (np.zeros( (Height,Width*(Width-1)) ) 
                                                  for _ in range(8))
    
    dT = np.roll(t,-1,0) - t
    dT = dT[:-1].reshape(len(dT)-1,1)

    j = 0
    BLc, BLrc = np.zeros(((Height*Width*(Width-1)),3)),np.zeros(((Height*Width*(Width-1)),3))
    for i in range(0,Width-1):
        BLx[:,j:j+Width-i-1] = x[:,i].reshape((Height,1)) - x[:,(i+1)::]
        BLx[:, n:] = -BLx[:,:n] # redundant operants?
        BLy[:,j:j+Width-i-1] = y[:,i].reshape((Height,1)) - y[:,(i+1)::]
        BLy[:, n:] = -BLy[:,:n]
        BLz[:,j:j+Width-i-1] = z[:,i].reshape((Height,1)) - z[:,(i+1)::]
        BLz[:, n:] = -BLz[:,:n]
        
        blM[:,j:j+Width-i-1] = np.sqrt(BLx[:,j:j+Width-i-1]**2+BLy[:,j:j+Width-i-1]**2
                                       + BLz[:,j:j+Width-i-1]**2)
        blM[:,n:] = blM[:,:n] #blM[:,i:Width-i-1]
        BLrx[:,j:j+Width-i-1] = vx[:,i].reshape((Height,1)) - vx[:,(i+1)::]
        BLrx[:, n:] = -BLrx[:,:n]
        BLry[:,j:j+Width-i-1] = vy[:,i].reshape((Height,1)) - vy[:,(i+1)::]
        BLry[:, n:] = -BLry[:,:n]
        BLrz[:,j:j+Width-i-1] = vz[:,i].reshape((Height,1)) - vz[:,(i+1)::]
        BLrz[:, n:] = -BLrz[:,:n]

        
        blrM[:,j:j+Width-i-1] = np.sqrt(BLrx[:,j:j+Width-i-1]**2+BLry[:,j:j+Width-i-1]**2
                                        + BLrz[:,j:j+Width-i-1]**2)
        blrM[:,n:] = blrM[:,:n]#blrM[:,i:Width-i]
        j += Width-i-1
    
    BLc[:,0] = BLx.flatten()
    BLc[:,1] = BLy.flatten()
    BLc[:,2] = BLz.flatten()
    
    BLrc[:,0] = BLrx.flatten()
    BLrc[:,1] = BLry.flatten()
    BLrc[:,2] = BLrz.flatten()
    
    # BLc = np.linalg.norm(BLc,axis=1)
    
    return BLx, BLy, BLz, BLc, BLrx, BLry, BLrz, BLrc, blM, blrM
    
def GetBaseLinesNoV(x,y,z, t):
    """Get the baselines from a series of points spread in cartesian three-dimensional space\n
    Returns baseline coordinates, velocities and magnitudes spread and in combined arrays: \n
    BLx, BLy, BLz, BLc, BLrx, BLry, BLrz, BLrc, blM, blrM"""
    Width = x.shape[1]
    Height = x.shape[0]
    n = np.int32(Width*(Width-1)/2)
    
    BLx, BLy, BLz, blM = (np.zeros( (Height,Width*(Width-1)) ) 
                                                  for _ in range(4))
    
    dT = np.roll(t,-1,0) - t
    dT = dT[:-1].reshape(len(dT)-1,1)

    j = 0
    BLc = np.zeros(((Height*Width*(Width-1)),3))
    for i in range(0,Width-1):
        BLx[:,j:j+Width-i-1] = x[:,i].reshape((Height,1)) - x[:,(i+1)::]
        BLx[:, n:] = -BLx[:,:n] # redundant operants?
        BLy[:,j:j+Width-i-1] = y[:,i].reshape((Height,1)) - y[:,(i+1)::]
        BLy[:, n:] = -BLy[:,:n]
        BLz[:,j:j+Width-i-1] = z[:,i].reshape((Height,1)) - z[:,(i+1)::]
        BLz[:, n:] = -BLz[:,:n]
        
        blM[:,j:j+Width-i-1] = np.sqrt(BLx[:,j:j+Width-i-1]**2+BLy[:,j:j+Width-i-1]**2
                                       + BLz[:,j:j+Width-i-1]**2)
        blM[:,n:] = blM[:,:n] #blM[:,i:Width-i-1]
        
        j += Width-i-1
    
    BLc[:,0] = BLx.flatten()
    BLc[:,1] = BLy.flatten()
    BLc[:,2] = BLz.flatten()
    

    
    # BLc = np.linalg.norm(BLc,axis=1)
    
    return BLx, BLy, BLz, BLc,  blM, 

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

def sampleSphereGrid(aStep = 1.0, lonStart = 0.0, lonStop = 180.0 , latStart = 0.0, latStop = 180.0):
    """Create a set of points along a sphere, distributed following a grid lattice \n
    Returns a set of 3-dimensional cartesian points."""

    
    longitude = np.arange(lonStart, lonStop+aStep, aStep)
    stepslon = len(longitude)
    longitude = longitude.reshape(1,stepslon)
    longitude *= np.pi/180
    
    latitude = np.arange(latStart, latStop+aStep, aStep)
    stepslat = len(latitude)
    latitude = latitude.reshape(1,stepslat)
    latitude *= np.pi/180
    
    N = stepslon * stepslat
    
    normals = np.zeros((3,N))
    
    
    for i in range(0,latitude.shape[1]):
        
        normals[0, i*stepslon:(i+1)*stepslon] =  np.cos(longitude)*np.cos(latitude[0,i]) 
        normals[1, i*stepslon:(i+1)*stepslon] =  np.sin(longitude)*np.cos(latitude[0,i])
        normals[2, i*stepslon:(i+1)*stepslon] =  np.sin(latitude[0,i])
    
    return normals

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

def NormalizeArrayToDtype(arr, dType):
    
    warn.filterwarnings("ignore")
    im = np.array(arr/np.max(arr) * 255 , dtype = dType)
    return im

#%% PSF manipulation and cost functions
def costFunction(psf):
    """Simple cost function for coarse PSF grids"""
    W,H = psf.shape
    
    weights = np.zeros(psf.shape)
    
    # x,y = np.mgrid[-2:2:W*j, -2:2:H*j]
    
    # Radial distance weights
    # offset = 2/300
    # x,y = np.mgrid[-300:300-offset:300j, -300:300-offset:300j]
    # weights = np.sqrt(x**2+y**2)
    
    # Square box weights
    weights = np.ones(psf.shape)
    
    weights[int(W/2),int(H/2)] = 0
    cost = np.sum(psf*weights)/(W*H)
    # cost = (np.sum(psf) -psf[int(W/2),int(H/2)])/(W*H)
    
    return cost

from scipy.special import j1 as bessel
from os import getcwd
from PlottingFunctions import Simple2DImage
def costFunctionAiry(psf):
    """Simple cost function for coarse PSF grids"""
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
    # cost = (np.sum(psf) -psf[int(W/2),int(H/2)])/(W*H)
    
    # cwd = getcwd()
    # Simple2DImage(NormalizeArrayToDtype(np.abs(np.fft.fftshift(airy)), np.uint8),
    #               title=f"Airy function reference", savefig = True, 
    #               figfolder = cwd, name = 'Airy.png', xlabel = 'x [km]', ylabel = 'y [km]', colormap = 'Greys_r',
    #               axes = 'psf')
    return cost

def PSFanalysis(normals, samplepoints, resolution = 100., cutout = 0, size=100e3):
    """Do a analysis of PSF for a series of view directions represented by a array of normal vectors \n
    Returns a complex128 poins spread function, and the uv sampling grid"""
    # 
    
    L = int(size/resolution)    
    N = normals.shape[1]
    cost = np.zeros(N)
    
    grid, I_psf = np.zeros((2*L,2*L,N)),np.zeros((2*L,2*L,N),dtype=np.complex128)
    
    for i in range(0,N):
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
 
        # # Clear up as much memory as possible for the next step
        # del v,u, Projection
        
        I_psf[:,:,i] = np.fft.ifft2(grid[:,:,i])
        
        cost[i] = costFunctionAiry(NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_psf[:,:,i])), np.uint8))
        
        # cwd = getcwd()
        # Simple2DImage(NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_psf[:,:,i])), np.uint8),
        #           title=f"Point spread function, cost: {cost[i]}", savefig = True, 
        #           figfolder = cwd, name = f"cost_{i}_10.png", xlabel = 'x [km]', ylabel = 'y [km]', colormap = 'Greys_r',
        #           axes = 'psf')
        # Simple2DImage(NormalizeArrayToDtype(np.abs(grid[:,:,i]), np.uint8),
        #           title="Sample function", savefig = True, 
        #           figfolder = cwd, name = f"sample_{i}_10.png", colormap = 'Greys_r',
        #           axes = 'baselines')
    
    if cutout!=0:
        Lc = int(cutout/resolution)    
        I_psf = I_psf[(L-Lc):(L+Lc), (L-Lc):(L+Lc)]
    
    return I_psf,grid,cost
        
    
def PSFanalysisSlow(normals, samplepoints,  resolution = 15., cutout = 0, size=100e3, printprogress = False):
    """Stepwise PSF analysis for a series of view directions represented by a array of normal vectors \n
    Returns a complex128 poins spread function, and the uv sampling grid"""
    # Do a analysis of PSF for a series of view directions represented by a array of normal vectors
    # NOT SEEING ANY ADVANTAGES OVER THE NORMAL IFFT2 ATM....
    L = int(size/resolution)    
    N = normals.shape[1]
    
    for i in range(0,N):
        Projection = Project2Dplane(normals[:,i], samplepoints)
        if printprogress:
            print("Now at normal index " + str(i) + " of " + str(N))

        x = Projection[0,:]/(resolution) + L # Was offset by +L to centralize the uv plane axis
        y = Projection[1,:]/(resolution) + L
        x = x.astype(int)
        y = y.astype(int)
        
        x[np.abs(y)>=2*L] = 3*L
        y[np.abs(y)>=2*L] = 3*L
        y[np.abs(x)>=2*L] = 3*L
        x[np.abs(x)>=2*L] = 3*L
        
        holder = np.zeros((2*L,2*L), dtype=np.complex64)
        grid = np.zeros((2*L,2*L), dtype=np.uint8)
        
        for rw in range(0,2*L):
            row = np.zeros((2*L,1))
            row[x[y==rw],0] = 1
            grid[rw,:] = row.flatten()
            holder[rw,:] = np.fft.ifft(row.flatten())
            
        for cl in range(0,2*L):
            holder[:,cl] = np.fft.ifft(holder[:,cl]).flatten()

    I_psf = holder 
    
    if cutout!=0:
        Cutrange = cutout
        Lc = int(Cutrange/resolution)    
        I_psf = I_psf[(L-Lc):(L+Lc), (L-Lc):(L+Lc)]
    
    grid = grid.astype(np.uint8)
    return I_psf, grid

def ImagingAnalysis(uv_dist, I_sky,figfolder,resolution = 1, figadd = ''):
    """Do a full imaging analysis of a reference sky image from a baseline set""" 
    I_sky += + 1 *(np.abs(I_sky) == 0).astype(np.uint8) 
    
    from matplotlib import pyplot as pp
    
    import PlottingFunctions as PF
    
    fig = PF.Simple2DImage(I_sky,title="Original sky image", xlabel = "pixels", ylabel = "pixels")
    PF.saveFig(figfolder,'I_sky'+ figadd + '.png')
    
    
    # fig = pp.figure(figsize=(4,4))
    # pp.title('Original sky image')
    # pp.axis([-150,150,-150,150])
    # pp.imshow(I_sky,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    # pp.savefig(figfolder + 'I_sky'+ figadd + '.png')
    
    L = np.floor(I_sky.shape[0]/2)
    I_fft = np.fft.fft2(I_sky)
    
    fft_shift = np.fft.fftshift(I_fft)
    
    im_fft = np.abs(fft_shift)
    im_fft = NormalizeArrayToDtype(im_fft,np.uint8)

    # fig = pp.figure(figsize=(4,4))
    # pp.title('Centralized fourier transform of sky image')
    # pp.axis([-150,150,-150,150])
    # pp.imshow(im_fft,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    # pp.savefig(figfolder + 'CentralFourier'+ figadd +'.png')

    fig = PF.Simple2DImage(im_fft,title="Centralized fourier transform of sky image",
                              axes = "psf")
    PF.saveFig(figfolder,'CentralFourier'+ figadd +'.png')
    
    # size = np.int(np.nanmax(np.abs(uv_dist)/(resolution) + L))
    
    V = np.zeros(I_sky.shape,dtype=np.complex128)
    # V = np.zeros((size,size), dtype = np.complex64)
    for i in range(0,uv_dist.shape[1]):
        u = np.int(uv_dist[0,i]/(resolution) + L) 
        v = np.int(uv_dist[1,i]/(resolution) + L)
        
        # Why were these flipped?
        # V[v,u] = fft_shift[u,v]
        V[u,v] = fft_shift[u,v]
        #V[u,v] = I_fft[u,v]
        
    
    S = np.zeros(V.shape) + 1.* (np.abs(V) > 0) 
    I_PSF = NormalizeArrayToDtype(np.abs(
        np.fft.ifftshift(np.fft.ifft2(
            (S.astype(np.complex128))
                      ))),np.uint8)
    
    # fig = pp.figure(figsize=(4,4))
    # pp.title('Point Spread Function')
    # pp.axis([-150,150,-150,150])
    # pp.imshow(I_PSF,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    # pp.savefig(figfolder + 'PSF'+ figadd +'.png')
    
    fig = PF.Simple2DImage(I_PSF,title="Point Spread Function",
                              axes = "psf")
    PF.saveFig(figfolder,'PSF'+ figadd +'.png')
    
    I_im = np.fft.ifft2(np.fft.ifftshift(V))
    
    # fig = pp.figure(figsize=(4,4))
    # pp.title('Sample function')
    # pp.axis([-150,150,-150,150])
    # pp.imshow(S,origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    # pp.savefig(figfolder + 'SampleFunc'+ figadd +'.png')
    
    fig = PF.Simple2DImage(S,title="Sample function",
                              axes = "baselines")
    PF.saveFig(figfolder,'SampleFunc'+ figadd +'.png')
    
    # fig = pp.figure(figsize=(4,4))
    # pp.title('"Dirty" sky image')
    # pp.axis([-150,150,-150,150])
    # pp.imshow(NormalizeArrayToDtype(I_im, np.uint8),origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    # pp.savefig(figfolder + 'I_obs'+ figadd +'.png')
    
    fig = PF.Simple2DImage(NormalizeArrayToDtype(np.abs(I_im), np.uint8),title='"Dirty" sky image',
                              xlabel = "pixels", ylabel = "pixels")
    PF.saveFig(figfolder,'I_obs'+ figadd +'.png')
    
    # fig = pp.figure(figsize=(4,4))
    # pp.title('Centralized fourier of dirty image')
    # pp.axis([-150,150,-150,150])
    # pp.imshow(NormalizeArrayToDtype(np.fft.fftshift(np.fft.fft2(I_im)), np.uint8),origin="lower", cmap='Greys_r', extent=(-150,150,-150,150))
    # pp.savefig(figfolder + 'I_fourier'+ figadd +'.png')
    
    fig = PF.Simple2DImage(NormalizeArrayToDtype(np.fft.fftshift(np.fft.fft2(I_im)), np.uint8),
                              title='Centralized fourier of dirty image',
                              axes="psf")
    PF.saveFig(figfolder,'I_fourier'+ figadd +'.png')
    
    return im_fft, V, I_im, I_PSF

def ConvertL4Coordinates(pointstate, centerstate = []):
    """Compute the L4 point location in a cartesian reference frame around a central body \n
    Returns a set of cartesian coordinates for the L4 point"""
    
    P_m = pointstate[:,0:3]
    V_m = pointstate[:,3:]
    
    if V_m == []:
        print("Missing velocity data for L4 conversion!")
        return []
    
    if centerstate == []:
        centerstate = np.zeros(P_m.shape)
    else:
        centerstate = centerstate[:,:3]
    
    P_m -= centerstate
    
    N = len(pointstate[:,1])
    H_m = np.cross(P_m, V_m)
    L4dir = -np.cross(P_m , H_m)
    
    L4 = centerstate + 0.5*P_m + (0.5*np.sqrt(3) *  np.linalg.norm(P_m,axis=1).reshape(N,1)
        * np.ones((N,3)) * L4dir / np.linalg.norm(L4dir,axis=1).reshape(N,1) * np.ones((N,3)))

    return L4

def ExtrapolateSwarmMotion(x, y, z, pos):
    N = x.shape[0]
    # Propagate motion of the core:
    
    delta_p = np.array( (np.mean(np.roll(x,-1,axis=0)-x,axis=1),
                         np.mean(np.roll(y,-1,axis=0)-y,axis=1),
                         np.mean(np.roll(z,-1,axis=0)-z,axis=1) ))
    delta_p[:,-1] = np.zeros((3))
    delta_p = np.roll(delta_p.T, 1, axis=0)
    
    extra_pos = np.zeros((N,3))
    extra_pos[0,:] = pos
    
    for i in range(1,N):
        extra_pos[i,:] = extra_pos[i-1,:] + delta_p[i,:]
        
    return extra_pos

def vectorProj(v1,v2):
    """Projection of vector v1 onto v2"""
    
    v3 = np.dot(v1,v2)/np.linalg.norm(v2)
    return v3
