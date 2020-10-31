# -*- coding: utf-8 -*-
"""
Cost function creation / test suite
Created on Thu May  7 11:42:27 2020

@author: USER
"""


import numpy as np
import sys 

from os import chdir
from os import getcwd


sys.path.append('../')
chdir('../')
import LazyLib as LL
import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF

figure_folder = LL.folder_dict_funcs["CF"]


print("Start synchronizing data files")
# LL.SyncDataFiles()
print("Synchronized data files")
cwd = getcwd()
datafolder = cwd + '/Data/PerturbationAnalysis/'
File = datafolder + 'propagationHistory_0_fullprop.dat'

print("Importing data")
t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,1, True)
print("Imported data succesfully")
figfolder = cwd + '/Figures/Cost Function/'
LL.CheckDir(figfolder)

chdir('Cost function')

PF.CloseAll()
#%%

sat = np.concatenate((x,y,z), axis=1)
co_x , co_y, co_z, L4b = DI.ConvertToBarycentricFrame(moon, sat, exportL4 = True, fixedL4 = True)
baryfig, baryax = PF.Simple3DPlot(co_x,co_y,co_z)
PF.Simple3DScatter(L4b[-1,0],L4b[-1,1],L4b[-1,2], fig=baryfig, ax=baryax)
PF.Simple3DPlot(L4b[:,0], L4b[:,1], L4b[:,2], 'Barycentric frame motion', xlim=(-0,0.4E6), ylim=(-0,0.6E6), fig=baryfig, ax = baryax)

N = co_x.shape[0]
fig2, ax2 = PF.Simple3DPlot(co_x - L4b[:,0].reshape(N,1)*np.ones(co_x.shape), co_y - L4b[:,1].reshape(N,1)*np.ones(co_y.shape)
                , co_z - L4b[:,2].reshape(N,1)*np.ones(co_z.shape))
PF.Simple3DScatter(0,0,0, fig= fig2, ax = ax2, title="Relative motion to L4 in barycentric frame")

#%% 
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = False, baselineMagnitude = BL_m, 
                   baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)
PF.saveFig(figure_folder, "3D baselines")

L4 = GF.ConvertL4Coordinates(moon)
l4_x, l4_y, l4_z = DI.ConvertToRelativeFrame(L4[:,0], L4[:,1], L4[:,2], x,y,z)

PF.Simple3DPlot(l4_x, l4_y, l4_z, title="L4-centralized satellite positions")
#%% 
# normals = GF.sampleSphereGrid(aStep = 6, latStop = 180, lonStop = 179)
# PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], title='Grid distributed sample points along unit sphere')
normals = GF.sampleSphereFib(30)
PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], 'Fibonnaci lattice sample points along unit sphere')

#%%

# Discarding invalid baselines is already done internally!
# BLm = np.linalg.norm(BL_c, axis=1)
# for i in range(0,BL_c.shape[0]-1):
#     if BLm[i] >= 100e3:
#         BL_c[i,:] = np.nan*np.ones((1,3))


res = 200e3/50

I_PSF, I_S,cost = GF.PSFanalysis(normals, BL_c,resolution = res)

index = int(np.where(cost == np.min(cost))[0])

WL = 30

#%%
PF.Simple2DImage(I_S[:,:,index], title='Sample function',axes = 'baselines')
PF.Simple2DImage(GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF[:,:,index])), np.uint8),
                 title='Point Spread Function',axes='psf')

#%% do the porcupine plot!
fig, ax = PF.BaselinePlot3DWithCost(BL_x,BL_y,BL_z,cost,normals, savefig = False, figfolder = '', name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                 scalemod=1e3)


# PF.Animate3DPlotToGIF(fig, ax, "Animated cost porcupine", figfolder)
#%% Experimenting with the wonky ass 3-d plot drawer Can always revert to scatter plot

cost_mean = np.mean(cost)


#%% Compare interior results of different resolutions

normals = np.array([[1],
                    [0],
                    [0]])


pixel_array = np.array([6800,3400,1200, 800, 600, 400, 200, 100, 50,10])
res0 = 200e3/6800
I_PSF50, I_S50,cost50 = GF.PSFanalysis(normals, BL_c,resolution = res0)

#%% Comparison of different resolutions and their effects on the PSF
cw = 25
angle_width = cw*np.degrees(np.sin(1.22*30/100e3))*60
rms, sim,res = np.zeros(len(pixel_array)),np.zeros(len(pixel_array)),np.zeros(len(pixel_array))
pix_red = np.empty_like(rms)
deviation = np.zeros((2*cw,2*cw,len(pixel_array)))

rms[0] = 0
res[0] = res0
sim[0] = 100
pix_red[0] = 0

acorr = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSF50)), np.int)
Nr = int(I_PSF50.shape[0]/2)
PF.Simple2DImage(acorr[Nr-cw:Nr+cw,Nr-cw:Nr+cw,0],
                      title=f"PSF, {round(res[0])}x{round(res[0])}m pixels",
                      savefig = True, figfolder = figure_folder,name = f"PSF_{pixel_array[0]}_cw",axes = "psf")

for i in range(1,len(pixel_array)):
    res[i] = 200e3/pixel_array[i]
    
    pix_red[i] = (1-(pixel_array[i]**2/6800**2))*100
    
    I_PSFref, I_Sref,costref = GF.PSFanalysis(normals, BL_c,resolution = res[i])
    
    
    bcorr = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFref)), np.int) 
    Nc = int(pixel_array[i]/2)
    
    error = acorr[Nr-cw:Nr+cw,Nr-cw:Nr+cw,0] - bcorr[Nc-cw:Nc+cw,Nc-cw:Nc+cw,0]
    
    rms[i] = np.sum(np.abs(error))
    
    PF.Simple2DImage(bcorr[Nc-cw:Nc+cw,Nc-cw:Nc+cw,0],
                      title=f"PSF, {round(res[i])}x{round(res[i])}m pixels",
                      savefig = True, figfolder = figure_folder,name = f"PSF_{pixel_array[i]}_cw",axes = "psf")
    
    PF.Simple2DImage(np.abs(error).astype(np.uint8),
                      title=f"PSF error, {round(res[i])}x{round(res[i])}m pixels",
                      savefig = True, figfolder = figure_folder,name = f"PSF_error_{pixel_array[i]}_cw", axes="psf")
    
    deviation[:,:,i] = np.abs(error)/255
    sim[i] = (1-rms[i]/(cw**2*255))*100

#%%
# Remove too large elements from bl array
import cv2
BLm = np.linalg.norm(BL_c, axis=1)
for i in range(0,BL_c.shape[0]-1):
    if BLm[i] >= 100e3:
        BL_c[i,:] = np.nan*np.ones((1,3))

BL_filtered = BL_c[np.isnan(BL_c) == False]
BL_filtered = BL_filtered.reshape(int(len(BL_filtered)/3),3)

l1 = 666.66*np.arange(-150,150)

spoof = np.zeros((300**2,2))
for i in range(0,300):
    spoof[300*i:300*(i+1),0] = l1[i]
    spoof[300*i:300*(i+1),1] = l1

refimg = cv2.imread(figfolder + 'Refimg_Andromeda.jpg',cv2.IMREAD_GRAYSCALE)
xvec = np.array([[1],[0],[0]])
projection= GF.Project2Dplane(xvec[:,0], BL_filtered)
im_fft, V, I_im, I_PSF = GF.ImagingAnalysis(projection, refimg,figfolder,666.66, 'andromeda_x')
# I_PSFx2, dump = GF.PSFanalysisSlow(xvec, BL_c,  resolution = 50.)
# I_PSFx2 = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFx2)), np.uint8)
# PF.Simple2DImage(I_PSFx2[1900:2100,1900:2100], title='PSF evaluated from the X-axis',axes = 'psf')

#%% Airy disk graph

from scipy.special import j1 as bessel
from matplotlib import pyplot as pp

figfolder = LL.folder_dict["Methodology"]

W,H = 500,500
ang = 3*np.degrees(np.arcsin(1.22*30/100e3))/500 #deg / px
angles = np.zeros((W,H))
x,y = np.arange(-W/2*ang,W/2*ang,ang),np.arange(-H/2*ang,H/2*ang,ang)
x = x.reshape(1,len(x))
y = y.reshape(len(y),1)

WL = LL.Constants["c"]/10e6
R = 50e3

angles = np.sqrt(x**2 + y**2)
x = 2*np.pi/WL * R * np.sin(np.deg2rad(angles))
airy = (255*( 2* bessel(x)/x)**2).astype(int)

exw = W/2*ang * 60
exh = H/2*ang * 60

fig = pp.figure(figsize=(5,5))
pp.title("Airy function at 10 MHz, B=100km")
im = pp.imshow( GF.NormalizeArrayToDtype(airy,np.uint8), origin="lower", 
                  cmap="pink", extent=(-exw,exw,-exh, exh))
pp.xlabel('u [arcsmin]')
pp.ylabel('v [arcmin]')
pp.colorbar(im,fraction=0.046, pad=0.04)

pp.tight_layout()
PF.saveFig(figfolder, "Airy function")