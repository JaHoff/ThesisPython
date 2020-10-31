# -*- coding: utf-8 -*-
"""
Container for the functions used in the verification algorithm
Created on Wed Mar 11 10:34:16 2020

@author: Jurriaan van ' Hoff
"""

from os import getcwd
from os import chdir
import sys 
import numpy as np
from matplotlib import pyplot as pp 

sys.path.append('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF


#%% DATA IMPORT FUNCTIONS
def VerifyRelativeFrame():

    f_x = np.array([[0], [1]])
    f_y = np.zeros((2,1))
    f_z = np.zeros((2,1))
    
    x = np.array([[3.], [3.]])
    y = np.array([[3.], [3.]])
    z = np.array([[3.], [3.]])
    
    d = LL.RelativeFolder("..//Figures/Verification/")
    fig,ax = PF.PlotAnnotated3D(f_x.flatten(), f_y.flatten(), f_z.flatten(), lgnd = 'Frame origins')
    PF.AddVector(ax,[f_x[1],f_y[1],f_z[1]], [x[1],y[1],z[1]])
    PF.PlotAnnotated3D(x.flatten(), y.flatten(), z.flatten(), title='Relative frame test case',
                       savefig = True, figfolder = d, name = 'RelFrame', fig=fig, ax=ax, mkr='x', lgnd = 'Point coordinates')
    pp.show()
    
    
    xf,yf,zf = DI.ConvertToRelativeFrame (f_x, f_y, f_z, x, y, z)
  
    
    
    # expected output
    xe = np.array([[3.], [2.]])
    ye = np.array([[3.], [3.]])
    ze = np.array([[3.], [3.]])

    if np.all(xe == xf) and np.all(ye == yf) and np.all(ze == zf):
        return True
        
    
    return False

def VerifyBarycentric():
    
    body2_p = np.array([ [10., 0., 0., 0., 10., 0.],    
                         [10., 10., 0., 0., 10., 0.]])
    body3_p = np.array([ [5, 5., 5.],
                         [5., 5., 8.]])
    
    xf,yf,zf,L4 = DI.ConvertToBarycentricFrame(body2_p, body3_p, body2_m = 1, body1_p = None, body1_m = 9, exportL4 = True)
    
    ans = np.concatenate((xf,yf,zf), axis=1)
    anse = np.array([ [4., 5., 5.],
                      [np.sqrt(50) - 0.1*np.sqrt(200), 0., 8.]])
    
    ans=np.round(ans,6)
    anse = np.round(anse,6)
    
    
    body1_p = np.array([ [0., 0., 0., 5., 5., 0.],
                         [5., 5., 0., 5., 5., 0.]])
    body2_p = np.array([ [10., 0., 0., 0., 10., 0.],
                         [10., 10., 0., 0., 10., 0.]])
    body3_p = np.array([ [5, 5., 5.],
                         [5., 5., 8.]])
    
    xf,yf,zf,L4 = DI.ConvertToBarycentricFrame(body2_p, body3_p, body2_m = 1, body1_p = body1_p, body1_m = 9, exportL4 = True)
    
    ans2 = np.concatenate((xf,yf,zf), axis=1)
    anse2 = np.array([ [4., 5., 5.],
                      [- 0.1*np.sqrt(50), 0., 8.]])
    
    # to avoid failure on numerical residuals
    ans2=np.round(ans2,6)
    anse2 = np.round(anse2,6)
    
    if np.all((ans == anse)) and np.all((ans2 == anse2)):
        return True
    
    
    return False

#%% GENERAL FUNCTIONS
    
def VerifyBaselines():
    
    t = np.array([[0],
                  [1]])
    
    x = np.array([[1 , 3],
                  [-4, 3]])
    y = np.array([[0 , 1],
                  [6, -5]])
    z = np.array([[0 , 0],
                  [-4, 3]])
                 
    vx = np.array([[-5 , 0],
                   [-5, 0]])
    vy = np.array([[6 , -6],
                   [6, -6]])
    vz = np.array([[-4 , 3],
                   [-4, 3]])
    
    BLx, BLy, BLz, BLc, BLrx, BLry, BLrz, BLrc, blM, blrM = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
    
    BLxa = np.array([[-2,2],
                     [-7,7]])
    BLya = np.array([[-1,1],
                     [11,-11]])
    BLza = np.array([[-0,0],
                     [-7,7]])
    BLrxa = np.array([[-5,5],
                     [-5,5]])
    BLrya = np.array([[12,-12],
                     [12,-12]])
    BLrza = np.array([[-7,7],
                     [-7,7]])
    
    M0 = np.sqrt(BLxa[0,0]**2 + BLya[0,0]**2 + BLza[0,0]**2)
    M1 = np.sqrt(BLxa[1,0]**2 + BLya[1,0]**2 + BLza[1,0]**2)
    blMa = np.array([[M0,M0],
                     [M1,M1]])
    M0 = np.sqrt(BLrxa[0,0]**2 + BLrya[0,0]**2 + BLrza[0,0]**2)
    M1 = np.sqrt(BLrxa[1,0]**2 + BLrya[1,0]**2 + BLrza[1,0]**2)
    blrMa = np.array([[M0,M0],
                     [M1,M1]])
    
    if( np.all(BLx == BLxa) and np.all(BLy == BLya)  and np.all(BLz == BLza)  and np.all(BLrx == BLrxa) 
       and np.all(BLry == BLrya)  and np.all(BLrz == BLrza)  and np.all(blM == blMa)  and np.all(blrM == blrMa) ):
        return True
        
    
    
    return False

def Verify2DProj():

    data = np.array([[1,1,1],
                     [0, 6, -7],
                     [3,-4, -5]])
    
    normal = np.array([0,7,0])   
    projection1 = GF.Project2Dplane(normal, data)
    
    ans1 = np.array([[ 1., -7., -5.],
                     [ 1.,  0.,  3.]
                     ])
    
    normal = np.array([1,0,0])   
    projection2 = GF.Project2Dplane(normal, data)
    
    ans2 = np.array([[ 1., -7., -5.],
                     [ -1.,  -6.,  4.]
                     ])
    
    normal = np.array([0,0,5])   
    projection3 = GF.Project2Dplane(normal, data)
    
    ans3 = np.array([[ 1., 0., 3.],
                     [ 1.,  6.,  -4.]])

    if np.all(projection1 == ans1) and np.all(projection2 == ans2) and np.all(projection3 == ans3):
        return True
    
    return False

def VerifyPSF():
    d = LL.RelativeFolder("..//Figures/Verification/")
    
    normals = np.array([1,0,0]).reshape(3,1)
    samplepoints = np.array([[-2,-1,0],
                             [2, 1, 0]])
    res = 0.1
    sz = 10
    I_psf,grid,ccs = GF.PSFanalysis(normals, samplepoints, resolution = res, cutout = 0,size=sz)
    
    # PF.Simple2DImage(grid, title='Sample function',savefig = True, figfolder = d, 
    #                  name='SampleFunction', xlabel='', ylabel='')
    # PF.Simple2DImage(GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_psf)),np.uint8), 
    #                  title='Point Spread Function',savefig = True, figfolder = d,
    #                  name='ReconstructedPSF', xlabel='', ylabel='')
    
    
    sp = np.fft.fft(I_psf[:,0,0]).T
    freq = np.fft.fftshift(np.fft.fftfreq(200))
    freq = np.abs(freq[sp == np.max(sp)])
    
    if freq != (1/res)/(2*sz/res):
        print("Frequency reconstruction mismatch!")
        return False
    
    tester = np.zeros(sp.shape)
    tester[sp >= 1e-10] = 1
    if np.all(tester.astype(int) == grid[:,100,0]) == False:
        print("Indirect reconstruction mismatch")
        return False
    
    freq = np.fft.fftshift(np.fft.fftfreq(200))
    # pp.figure(figsize=(4,4))
    # pp.plot(freq, sp.real, freq, sp.imag)
    # pp.title('Frequency spectrum of horizontal slice')
    # pp.xlabel('Frequency [Hz]')
    # pp.ylabel('Magnitude')
    # pp.xlim((-0.1,0.1))
    # pp.savefig(d+'frequency spectrum.png')
    
    normals = np.random.rand(3,1)
    samplepoints = 10*np.random.rand(3,3)
    
    I_psf,grid,cc = GF.PSFanalysis(normals, samplepoints, resolution = res, cutout = 0,size=sz)
    
    sp = np.fft.fft2(I_psf[:,:,0])
    recongrid = GF.NormalizeArrayToDtype((sp),np.uint8)
    recongrid[recongrid > 0] = 255
    
    xg,yg = np.where(grid[:,:,0] > 0)
    
    xrg,yrg = np.where(recongrid > 0)
    if np.all(recongrid == GF.NormalizeArrayToDtype(grid[:,:,0],np.uint8)) == False:
        print("Complex grid reconstruction mismatch!")
        return False

    return True

def VerifyPSFSlow():
    normals = np.array([1,0,0]).reshape(3,1)
    samplepoints = np.array([[-2,-1,0],
                             [2, 1, 0]])
    res = 0.1
    sz = 10
    I_psf,grid = GF.PSFanalysisSlow(normals, samplepoints, resolution = res, cutout = 0,size=sz)
    
    PF.ClearPlots()
    PF.Simple2DImage(grid, title='Sample function')
    PF.Simple2DImage(GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_psf.T)),np.uint8), title='Point Spread Function')
    
    
    sp = np.fft.fft(I_psf[0,:])
    freq = np.fft.fftshift(np.fft.fftfreq(200))
    freq = np.abs(freq[sp == np.max(sp)])
    
    if freq != (1/res)/(2*sz/res):
        print("Frequency reconstruction mismatch!")
        return False
    
    tester = np.zeros(sp.shape)
    tester[sp >= 1e-10] = 1
    if np.all(tester.astype(int) == grid[100,:]) == False:
        print("Indirect reconstruction mismatch")
        return False
    
    
    # pp.figure()
    # pp.plot(freq, sp.real, freq, sp.imag)
    
    normal = np.random.rand(1,3)
    samplepoints = 10*np.random.rand(3,3)
    
    res = 0.1
    sz = 10
    I_psf,grid = GF.PSFanalysisSlow(normals, samplepoints, resolution = res, cutout = 0,size=sz)
    
    sp = np.fft.fft2(I_psf)
    recongrid = GF.NormalizeArrayToDtype((sp),np.uint8)
    recongrid[recongrid > 0] = 255
    
    if np.all(recongrid == GF.NormalizeArrayToDtype(grid,np.uint8)) == False:
        print("Complex grid reconstruction mismatch!")
        return False
    
    return True

def VerifyL4Coordinates():
    
    pointstate = np.array([[0,2,0,1,0,0],
                           [1,2,1, 1,0,0]])
    centerstate = np.array([[0,0,0],
                           [1,0,0]])
    L4 = GF.ConvertL4Coordinates(pointstate, centerstate)
    
    ans = np.array([[np.sqrt(3),1,0],
                    [1+np.sqrt(15)*0.5,1,0.5]])
    
    if np.all(ans != L4):
        return False
    return True


