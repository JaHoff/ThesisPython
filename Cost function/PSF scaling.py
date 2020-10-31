# -*- coding: utf-8 -*-
"""
Cost function creation / test suite
Created on Thu May  7 11:42:27 2020

@author: USER
"""


import numpy as np
import sys 
# import cv2
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
File = datafolder + 'propagationHistory_0_30km.dat'

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
# baryfig, baryax = PF.Simple3DPlot(co_x,co_y,co_z)
# PF.Simple3DScatter(L4b[-1,0],L4b[-1,1],L4b[-1,2], fig=baryfig, ax=baryax)
# PF.Simple3DPlot(L4b[:,0], L4b[:,1], L4b[:,2], 'Barycentric frame motion', xlim=(-0,0.4E6), ylim=(-0,0.6E6), fig=baryfig, ax = baryax)

N = co_x.shape[0]
# fig2, ax2 = PF.Simple3DPlot(co_x - L4b[:,0].reshape(N,1)*np.ones(co_x.shape), co_y - L4b[:,1].reshape(N,1)*np.ones(co_y.shape)
                # , co_z - L4b[:,2].reshape(N,1)*np.ones(co_z.shape))
# PF.Simple3DScatter(0,0,0, fig= fig2, ax = ax2, title="Relative motion to L4 in barycentric frame")

#%% 
BL_x, BL_y, BL_z, BL_c,  BLrx, BLry, BLrz, BLrc, BL_m, BLr_m  = GF.GetBaseLines(x,y,z, vx, vy, vz, t)
PF.UVWBaselinePlot(BL_x, BL_y, BL_z, includeNegatives=True, highlightfaulty = True, baselineMagnitude = BL_m, 
                   baselineRateMagnitude = BLr_m, baselineRateMagnitudeThreshold = 1)

L4 = GF.ConvertL4Coordinates(moon)
l4_x, l4_y, l4_z = DI.ConvertToRelativeFrame(L4[:,0], L4[:,1], L4[:,2], x,y,z)

PF.Simple3DPlot(l4_x, l4_y, l4_z, title="L4-centralized satellite positions")
#%% 
# normals = GF.sampleSphereGrid(aStep = 6, latStop = 180, lonStop = 179)
# PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], title='Grid distributed sample points along unit sphere')
normals = GF.sampleSphereFib(30)
# PF.Simple3DScatter(normals[0,:], normals[1,:], normals[2,:], 'Fibonnaci lattice sample points along unit sphere')

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


pixel_array = np.array([6800,3400,1200, 800, 600, 400, 200, 100, 50])
res0 = 200e3/6800
I_PSF50, I_S50,cost50 = GF.PSFanalysis(normals, BL_c,resolution = res0)

#%%
cw = 15
angle_width = np.round(cw*np.degrees(np.sin(1.22*30/100e3))*60,2)
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
                     title=f"Central section of PSF from {round(res[0])}x{round(res[0])}m pixels",
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
                     title=f"Central section of PSF from {round(res[i])}x{round(res[i])}m pixels",
                     savefig = True, figfolder = figure_folder,name = f"PSF_{pixel_array[i]}_cw",axes = "psf")
    
    PF.Simple2DImage(GF.NormalizeArrayToDtype(np.abs(error), np.uint8),
                     title=f"PSF error section from {round(res[i])}x{round(res[i])}m pixels",
                     savefig = True, figfolder = figure_folder,name = f"PSF_error_{pixel_array[i]}_cw", axes="psf")
    
    deviation[:,:,i] = np.abs(error)/255
    sim[i] = (1-rms[i]/(cw**2*255))*100

#%%
from matplotlib import pyplot as pp
fig,ax = PF.MakeFig(figsize=(12,3.5))
ax.plot(res,sim,color='b')
ax.plot(res,pix_red,color='r')
ax.scatter(res,sim,color='b',marker='d')
ax.scatter(res,pix_red,color='r',marker='d')
PF.DressFig(fig, ax, title=f"Effect of downscaling on central {angle_width} x {angle_width} arcmin of PSF",
            xlabel = "Sample pixel resolution [m]", ylabel = "Similarity/reduction [%]",ylim=[50,101],
            legendlist=['PSF similarity in %', 'Pixel reduction in %'])
pp.tight_layout()
PF.saveFig(figure_folder, f"Reduction_similarity_{cw}.png")


# verts = normals*cost

# id_prior = normals[2,:]#+np.tan(normals[0,:]/normals[1,:])
# temp = id_prior.argsort()
# ranks = np.empty_like(temp)
# ranks[temp] = np.arange(len(id_prior))

# verts_sort = verts[:,temp]

# Nd2 = len(ranks)/2
# # import matplotlib.pyplot as pp

# x = verts[0,:]
# y = verts[1,:]
# z = verts[2,:]



# fig,ax = PF.MakeFig(dimensions='3D')
# # ax.plot_trisurf(verts[0,:], verts[1,:], verts[2,:])
# # ax.plot_trisurf(x.flatten(),y.flatten(),z.flatten())
# mask = verts[2,:] >= 0
# mask = ranks >= Nd2-6
# ax.plot_trisurf(x[mask],y[mask],z[mask],color='b')
# mask = verts[2,:] <= 0
# mask = ranks <= Nd2 + 6
# ax.plot_trisurf(x[mask],y[mask],z[mask],color = 'b')
# PF.Simple3DSurf(verts[:,0], verts[:,1], verts[:,2], xmod = 1., ymod = 1., zmod = 1.,
#                  fig = None, ax = None)
#%% 


## WILL LIKELY NOT SEE CONCRETE DIFFERENCES IN THE PSF CORE UNLESS RESOLUTION IS AROUND A METER!
## USING 100 METER NOW YIELDS LITTLE DIFFERENCES WITH 10 METER IN GENERAL PIXEL SHAPE

# xvec = np.array([[1],[0],[0]])
# I_PSFx, I_Sx = GF.PSFanalysis(xvec, BL_c, resolution = 100)
# I_PSFx = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFx)), np.uint8)
# PF.Simple2DImage(I_Sx, title='Sample function x axis',axes = '100kmc')
# PF.Simple2DImage(I_PSFx[950:1050,950:1050], title='PSF function x axis',axes = '100kmc')

# xvec = np.array([[1],[0],[0]])
# I_PSFx, I_Sx,cost_x = GF.PSFanalysis(xvec, BL_c, resolution = 333.33)
# I_PSFx = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFx)), np.uint8)
# PF.Simple2DImage(I_Sx, title='Sample function x axis',axes = '100kmc')
# PF.Simple2DImage(I_PSFx, title='PSF function x axis',axes = '100kmc')


#%%
# Remove too large elements from bl array

# BLm = np.linalg.norm(BL_c, axis=1)
# for i in range(0,BL_c.shape[0]-1):
#     if BLm[i] >= 100e3:
#         BL_c[i,:] = np.nan*np.ones((1,3))

# BL_filtered = BL_c[np.isnan(BL_c) == False]
# BL_filtered = BL_filtered.reshape(int(len(BL_filtered)/3),3)

# l1 = 666.66*np.arange(-150,150)

# spoof = np.zeros((300**2,2))
# for i in range(0,300):
#     spoof[300*i:300*(i+1),0] = l1[i]
#     spoof[300*i:300*(i+1),1] = l1

# refimg = cv2.imread(figfolder + 'Refimg_Sixstars.jpg',cv2.IMREAD_GRAYSCALE)
# xvec = np.array([[0],[0],[1]])
# projection= GF.Project2Dplane(xvec[:,0], BL_filtered)
# im_fft, V, I_im, I_PSF = GF.ImagingAnalysis(projection, refimg,figfolder,666.66, '3star_z')
# I_PSFx2, dump = GF.PSFanalysisSlow(xvec, BL_c,  resolution = 50.)
# I_PSFx2 = GF.NormalizeArrayToDtype(np.abs(np.fft.fftshift(I_PSFx2)), np.uint8)
# PF.Simple2DImage(I_PSFx2[1900:2100,1900:2100], title='PSF x axis high resolution',axes = '100kmc')

