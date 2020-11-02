# -*- coding: utf-8 -*-
"""
Contains several sets of code used to ease plotting standard figure formats.

Code later deprecated and replaced by class-based plotting, see Presentation prep.py

As a result code is not fully maintained anymore

Old code still relies on this library however, which is why it is kept around

Created on Thu Feb 20 16:17:04 2020

@author: Jurriaan van 't Hoff
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as pp
import matplotlib.colors as col
import LazyLib as LL
from mpl_toolkits.mplot3d import Axes3D



font = {'family' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)

# %%3 dimensional plots of 3-D data:
def Simple3DPlot(x, y, z, title='', savefig = False, figfolder = '', name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                 xmod = 1000., ymod = 1000., zmod = 1000., xlim = None, ylim = None , zlim = None,
                 fig = None, ax = None):
    """Create a simple 3-d plot of a x,y,z dataset \n
    INPUT: \n
    x,y,z: Data arrays 1d or 2d, 2d will use [:,i] to plot multiple lines from a single column \n
    title: Title for figure \n
    savefig: Save figure? \n
    figfolder: folder to save figure in \n
    name: file name \n
    x/y/zlabel: labels to assign to x,y,z axes \n
    x/y/zmod: scale modifiers (def. 1000, converts to km) \n
    x/y/zlim: set limits for x,y,z axes \n
    fig,ax: Provide the function with an existing figure/axis ,otherwise create a new one.
    """
    
    if fig == None or ax == None:
        fig = pp.figure()
        ax = fig.add_subplot(111, projection='3d')
    
    if len(x.shape) == 2:  
        width = x.shape[1]
        for i in range(0,width):
            ax.plot(x[:,i]/xmod,y[:,i]/ymod,z[:,i]/zmod)
    else:
        ax.plot(x/xmod,y/ymod,z/zmod)
        
    # for i in range(0,width):
    #     ax.plot(x[:,i]/xmod,y[:,i]/ymod,z[:,i]/zmod)
    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)
    
    if savefig: 
        pp.savefig(figfolder + name)
    
#    fig.show()
    return fig, ax
    
def Simple3DScatter(x, y, z, title='', savefig = False, figfolder = '', name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]', 
                    zlabel='z [km]', xmod = 1000., ymod = 1000., zmod = 1000., fig = None, ax = None, color = None, marker=None):
    """Create a simple 3-d scatter plot of a x,y,z dataset \n
    INPUT: \n
    x,y,z: Data arrays 1d or 2d, 2d will use [:,i] to plot multiple lines from a single column \n
    title: Title for figure \n
    savefig: Save figure? \n
    figfolder: folder to save figure in \n
    name: file name \n
    x/y/zlabel: labels to assign to x,y,z axes \n
    x/y/zmod: scale modifiers (def. 1000, converts to km) \n
    x/y/zlim: set limits for x,y,z axes \n
    fig,ax: Provide the function with an existing figure/axis ,otherwise create a new one. \n
    color,marker: set scatterplot color and marker
    """
    if fig == None or ax == None:
        fig = pp.figure()
        ax = fig.add_subplot(111, projection='3d')
    try:
        if len(x.shape) == 2:  
            width = x.shape[1]
            for i in range(0,width):
                ax.scatter(x[:,i]/xmod,y[:,i]/ymod,z[:,i]/zmod,c=color, marker=marker)
        else:
            ax.scatter(x/xmod,y/ymod,z/zmod,c=color, marker=marker)
    except:
        ax.scatter(x/xmod,y/ymod,z/zmod,c=color, marker=marker)

    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    
    if savefig: 
        pp.savefig(figfolder + name)

    return fig, ax    


def PlotAnnotated3D(x,y,z,title='', savefig = False, figfolder = '', name = 'placeholder.png', xlabel = 'x', ylabel = 'y', zlabel='z',
                 xmod = 1., ymod = 1., zmod = 1., fig = None, ax = None, mkr=None, lgnd = None):
    """Create a simple 3-d scatter plot with annotations showing data value of a x,y,z dataset \n
    INPUT: \n
    x,y,z: Data arrays 1d or 2d, 2d will use [:,i] to plot multiple lines from a single column \n
    title: Title for figure \n
    savefig: Save figure? \n
    figfolder: folder to save figure in \n
    name: file name \n
    x/y/zlabel: labels to assign to x,y,z axes \n
    x/y/zmod: scale modifiers (def. 1000, converts to km) \n
    x/y/zlim: set limits for x,y,z axes \n
    fig,ax: Provide the function with an existing figure/axis ,otherwise create a new one. \n
    color,marker: set scatterplot color and marker
    """
    
    
    if fig == None or ax == None:
        fig = pp.figure()
        ax = fig.add_subplot(111, projection='3d')
    
    
    ax.scatter(x/xmod,y/ymod,z/zmod, marker = mkr,label=lgnd)
    
    for i in range(0,len(x)):
        label = '(%d, %d, %d)' % (x[i], y[i], z[i])
        ax.text(x[i], y[i], z[i], label)
        
    # if lgnd != None:
    #     ax.set_label(lgnd)

    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.legend()
    AddAxesWidget(ax)
    
    if savefig: 
        pp.savefig(figfolder + name)
    
    pp.show()
    return fig, ax 

def Simple3DSurf(x, y, z, xmod = 1000., ymod = 1000., zmod = 1000.,
                 fig = None, ax = None):
    """Create a simple 3-d surface plot of a x,y,z dataset \n
    INPUT: \n
    x,y,z: Data arrays [2d] \n
    x/y/zmod: scale modifiers (def. 1000, converts to km) \n
    fig,ax: Provide the function with an existing figure/axis ,otherwise create a new one. \n
    """
    if fig == None or ax == None:
        fig = pp.figure()
        ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x/xmod,y/ymod,z/zmod)
    
    return fig, ax

def Simple3DWireframe(x, y, z, xmod = 1000., ymod = 1000., zmod = 1000.,
                 fig = None, ax = None):
    """Create a simple 3-d wireframe plot of a x,y,z dataset \n
    INPUT: \n
    x,y,z: Data arrays [2d] \n
    x/y/zmod: scale modifiers (def. 1000, converts to km) \n
    fig,ax: Provide the function with an existing figure/axis ,otherwise create a new one. \n
    """
    if fig == None or ax == None:
        fig = pp.figure()
        ax = fig.add_subplot(111, projection='3d')

    ax.plot_wireframe(x/xmod,y/ymod,z/zmod)
    
    return fig, ax


#%% Two-dimensional data plots
def Simple2DImage(Image, title='', savefig = False, figfolder = '',
                  name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]',
                  colormap = 'Greys_r', axes = ''):
    """Create a simple 2-d image from a 2d array with color values \n
    INPUT: \n
    Image: image array \n
    title: Title for figure \n
    savefig: Save figure? \n
    figfolder: folder to save figure in \n
    name: file name \n
    x/ylabel: labels to assign to x,y axes \n
    colormap: colormap to use \n
    axes: Provide the function with an existing figure/axis ,otherwise create a new one. \n
    """
    
    L = Image.shape
    fig = pp.figure(figsize=(5,5))
    pp.title(title)
    
    
    if axes == '100kmc':
        pp.imshow( Image, origin="lower", cmap=colormap, extent=(-100,100,-100, 100))
        pp.xlabel('km')
        pp.ylabel('km')
    elif axes == 'baselines':
        ex = 100e3/30 / 1000
        pp.imshow( Image, origin="lower", cmap=colormap, extent=(-ex,ex,-ex, ex))
        pp.xlabel('k\u03BB u')
        pp.ylabel('k\u03BB v')
    elif axes == 'psf':
        
        W,H = Image.shape
        ang = np.degrees(np.arcsin(1.22*30/100e3))*60
        exw = W/2*ang
        exh = H/2*ang
        colormap = "pink"
        
        centers = [-exw,exw,-exh,exh]
        dx, = np.diff(centers[:2])/(H-1)
        dy, = -np.diff(centers[2:])/(W-1)
        extent = [centers[0]-dx/2, centers[1]+dx/2, centers[2]+dy/2, centers[3]-dy/2]

        pp.imshow( Image, origin="lower", norm = col.Normalize(vmin=0,vmax=255), 
                  cmap=colormap, extent=extent) #(-exw,exw,-exh, exh)
        step = exw/(int(W/10)*dx)
        pp.xticks(np.round(np.arange(centers[0], centers[1]+dx, step*dx),1)) #4.833333
        pp.yticks(np.round(np.arange(centers[3], centers[2]+dy, step*dy),1))
        pp.xlabel('u [arcmin]')
        pp.ylabel('v [arcmin]')
    else:
        pp.axis([-L[0],L[0],-L[1],L[1]])
        pp.imshow( Image, origin="lower", cmap=colormap, extent=(-L[0],L[0],-L[1],L[1]))
        pp.xlabel(xlabel)
        pp.ylabel(ylabel)
    
    
    if savefig: 
        pp.savefig(figfolder + name,  bbox_inches = "tight")

    return fig

def Simple2DPlot(x, y, title='', savefig = False, figfolder = '', name = 'placeholder.png',
                 xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]', xmod = 1000.,
                 ymod = 1000., xlim = None, ylim = None , fig = None, ax = None):
    """Create a simple 2-d plot from a x,y arrays \n
    INPUT: \n
    Image: image array \n
    title: Title for figure \n
    savefig: Save figure? \n
    figfolder: folder to save figure in \n
    name: file name \n
    x/ylabel: labels to assign to x,y axes \n
    colormap: colormap to use \n
    axes: Provide the function with an existing figure/axis ,otherwise create a new one. \n
    """
    if fig == None or ax == None:
        fig = pp.figure()
        ax = fig.add_subplot(11, projection='2d')
    
    if len(x.shape) == 2:  
        width = x.shape[1]
        for i in range(0,width):
            ax.plot(x[:,i]/xmod,y[:,i]/ymod)
    else:
        ax.plot(x/xmod,y/ymod)

    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if savefig: 
        pp.savefig(figfolder + name)
    
#    fig.show()
    return fig, ax

def Simple2DScatter(x, y, title='', savefig = False, figfolder = '', name = 'placeholder.png',
                 xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]', xmod = 1000.,
                 ymod = 1000., xlim = None, ylim = None , fig = None, ax = None, 
                 figsize = None, logScale = False,addAnnotation = False, 
                 annotation=None, color = None, legendAddition = None, texts = []):
    
    """Create a instant 2-dimensional scatter plot, or add scatter plots to an existing frame"""
    
    if fig == None or ax == None:
        fig,ax = pp.subplots(figsize=figsize)
    
    if len(x.shape) == 2:  
        width = x.shape[1]
        for i in range(0,width):
            ax.scatter(x[:,i]/xmod,y[:,i]/ymod, c = color, label=legendAddition)
    else:
        ax.scatter(x/xmod,y/ymod , c = color,label=legendAddition)

    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if addAnnotation:
        # if x < ax.get_xlim()[1] and x > ax.get_xlim()[0]:
        ax.text(x,y,annotation)
        texts.append(ax.text(x,y,annotation))
    
    if logScale:
        ax.set_xscale("log")
        ax.set_yscale("log")
    
    if legendAddition != None:
        ax.legend()
    
    if savefig: 
        pp.savefig(figfolder + name)
    
#    fig.show()
    return fig, ax

def Simple2DScatterAdjustedText(x, y, title='', savefig = False, figfolder = '', name = 'placeholder.png',
                 xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]', xmod = 1000.,
                 ymod = 1000., xlim = None, ylim = None , fig = None, ax = None, 
                 figsize = None, logScale = False,addAnnotation = False, 
                 annotation=None, color = None, legendAddition = None, texts = []):
    """Its like Simple2DScatter but the output of annotations is appended to a text array
    so that it is compatible with the adjustText package."""
    if fig == None or ax == None:
        fig,ax = pp.subplots(figsize=figsize)
    
    if len(x.shape) == 2:  
        width = x.shape[1]
        for i in range(0,width):
            ax.scatter(x[:,i]/xmod,y[:,i]/ymod, c = color, label=legendAddition)
    else:
        ax.scatter(x/xmod,y/ymod , c = color,label=legendAddition)

    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if addAnnotation:
        # if x < ax.get_xlim()[1] and x > ax.get_xlim()[0]:
        texts.append(ax.text(x,y,annotation))
    
    if logScale:
        ax.set_xscale("log")
        ax.set_yscale("log")
    
    if legendAddition != None:
        ax.legend()
    
    if savefig: 
        pp.savefig(figfolder + name)
    
#    fig.show()
    return fig, ax, texts
    
def Matrix2DPlot(x, M, title='', savefig = False, figfolder = '', name = 'placeholder.png',
                 xlabel = 'x [km]', ylabel = 'y [km]',  xmod = 1000., ymod = 1000.,
                 xlim = None, ylim = None , fig = None, ax = None, legendlist = None, 
                 stylelist=None, markN = 2,
                 plotAllEntries = False, logScale = False, figsize = None):
    """Plot the contents of a 3-dimensional matrix where M[x,0,0] yields an array 
    equal to time. \n plotAllEntries: Enable to plot all entries in axis 1, default only plots axis 0"""
    
    
    ## MOVE xmod scalingh ere! 
    sh = M.shape
    if fig == None or ax == None:
        fig = pp.figure(figsize=figsize)
        ax = fig.gca()#fig.add_subplot(11, projection= '2d')

    
    if (len(sh) != 3 and len(sh) !=2):
        print("Invalid matrix shape supplied!")
        
    

    if len(sh) == 2:
        for i in range(0,sh[1]):
            mask = np.isfinite(M[:,i])
            if stylelist != None:
                stl = stylelist[i]
                if stl != None:
                    ax.plot(x[mask]/xmod,M[mask,i]/ymod, stl, markevery = markN)
                    continue
            ax.plot(x[mask]/xmod,M[mask,i]/ymod)

    else:
        for i in range(0,sh[2]):
            if(plotAllEntries):
                for j in range(0,sh[1]):
                    mask = np.isfinite(M[:,j,i])
                    if stylelist != None:
                        stl = stylelist[i]
                        if stl != None:
                            ax.plot(x[mask]/xmod,M[mask,j,i]/ymod, stl, markevery = markN)
                            continue
                    ax.plot(x[mask]/xmod,M[mask,j,i]/ymod)
            else:
                mask = np.isfinite(M[:,0,i])
                if stylelist != None:
                        stl = stylelist[i]
                        if stl != None:
                            ax.plot(x[mask]/xmod,M[mask,0,i]/ymod, stl,  markevery = markN)
                            continue
                        
                ax.plot(x[mask]/xmod,M[mask,0,i]/ymod)
        
    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if logScale:
        ax.set_yscale('log')
        # ax.yaxis.set_major_locator(pp.MaxNLocator(20))
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    
    
    if legendlist != None:
        ax.legend(legendlist)
    
    if savefig: 
        pp.savefig(figfolder + name)
    
    return fig, ax


def Quiver2D(x,y,u,v,fig = None, color = None):
    
    pp.quiver(x,y,u,v,figure=fig, color = color)
    
    return

def UVWBaselinePlot(BLx,BLy,BLz, savefig = False, figfolder = '', name = 'placeholder.png', xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                 xmod = 1000., ymod = 1000., zmod = 1000., includeNegatives=True, highlightfaulty = False, baselineMagnitude = [],
                 baselineRateMagnitude = [],
                 baselineMagnitudeThreshold = 100e3, baselineRateMagnitudeThreshold = 5):
    fig = pp.figure()
    c = np.int(BLx.shape[1]/2)
    ax = fig.add_subplot(111, projection='3d')
    
    for i in range(0,c):
        ax.plot(BLx[:,i]/xmod, BLy[:,i]/ymod, BLz[:,i]/zmod, c='b', label = 'uvw baselines')
        
        if includeNegatives: 
            ax.plot(BLx[:,c+i]/xmod, BLy[:,c+i]/ymod, BLz[:,c+i]/zmod, c='g', label = 'negative pairs')
            
        if highlightfaulty:
            BLxplot, BLyplot, BLzplot = (np.zeros(BLx[:,i].shape) for _ in range(3))
            BLxplot[:] = BLx[:,i]
            BLxplot[baselineRateMagnitude[:,i]<baselineRateMagnitudeThreshold] = np.NaN
            BLyplot[:] = BLy[:,i]
            BLyplot[baselineRateMagnitude[:,i]<baselineRateMagnitudeThreshold] = np.NaN
            BLzplot[:] = BLz[:,i]
            BLzplot[baselineRateMagnitude[:,i]<baselineRateMagnitudeThreshold] = np.NaN
            ax.plot(BLxplot/xmod, BLyplot/ymod, BLzplot/zmod, c='y', linewidth=7.0)
            
            # BLxplot, BLyplot, BLzplot = [np.zeros(BLx[:,i].shape),]*3
            BLxplot[:] = BLx[:,i]
            BLxplot[baselineMagnitude[:,i]<baselineMagnitudeThreshold] = np.NaN
            BLyplot[:] = BLy[:,i]
            BLyplot[baselineMagnitude[:,i]<baselineMagnitudeThreshold] = np.NaN
            BLzplot[:] = BLz[:,i]
            BLzplot[baselineMagnitude[:,i]<baselineMagnitudeThreshold] = np.NaN
            ax.plot(BLxplot/xmod, BLyplot/ymod, BLzplot/zmod, c='r', linewidth=7.0)
            
            if includeNegatives:
                ax.plot(-BLxplot/xmod, -BLyplot/ymod, -BLzplot/zmod, c='r', linewidth=7.0)
            
    ax.set_title('uvw space baselines')    
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_zlabel('z [km]')
    if savefig:
        pp.savefig(figfolder + 'uvw3d.png')
        
    return fig,ax


def BaselinePlot3DWithCost(BLx,BLy,BLz,cost,normals, savefig = False, figfolder = '', name = 'placeholder.png', 
                           xlabel = 'x [km]', ylabel = 'y [km]', zlabel='z [km]',
                           scalemod=1e3):
    """Plot the baseline profile in 3D with a cost profile porcupine overlapped\n
    INPUT: \n
    BLx,BLy,BLz: [n,m] arrays with baseline data. n direction is time \n
    cost: [k] sized array with costs evaluated from a certain normal direction \n
    normals: [3,k] sized array with normal vectors"""
    
    fig,ax = MakeFig(figsize=(20,14), dimensions='3D')

    for i in range(0,BLx.shape[1]):
        ax.plot(BLx[:,i]/scalemod, BLy[:,i]/scalemod, BLz[:,i]/scalemod, c='b', label = 'uvw baselines')
            
    
    cost_mod = 1.1*np.max(np.abs([BLx,BLy,BLz]))/scalemod
    
    CostPorcupine3D(fig,ax,cost,normals,cost_mod)
    
    ax.set_title('uvw space baselines with evaluated cost function at sample points')    
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_zlabel('z [km]')
    if savefig:
        pp.savefig(figfolder + 'cost_porcupine.png')
        
    return fig,ax

#%% Plot cost function visualisers

def CostPorcupine3D(fig, ax, cost, normals,R):
    """Create a 3-dimensional porcupine plot using a set of direction vectors and magnitudes"""
    cost_corrected = (cost/np.mean(cost))**4
    verts = R*normals*cost_corrected
    
    
    ax.scatter(R*normals[0,:],R* normals[1,:], R*normals[2,:],c='k')
    # norm = col.Normalize(vmin=0.5,vmax=1.5)
    # ax.scatter(verts[0,:], verts[1,:], verts[2,:],c=cost_corrected, cmap = 'magma', norm=norm)
    
    for i in range(0,len(verts[0,:])):
        AddVector(ax,R*normals[:,i],verts[:,i], color='k')
    return

#%% QUALITY OF LIFE FUNCTIONS
def ClearPlots():
    """Close all figures"""
    pp.close("all")
    return

def CloseAll():
    """Close all figures\n
    Legacy function"""
    pp.close("all")
    return

def saveFig(folder, name):
    """Call savefig function outside of the normal constructor, in case this is
    more convenient"""
    LL.CheckDir(folder)
    pp.savefig(folder + name)
    return

def MakeFig(figsize = None, dimensions = '2D', proj = 'persp'):
    """Create a new figure, use dimensions to indicate 2D or 3D"""
    if dimensions == '2D':
        fig,ax = pp.subplots(figsize = figsize)
    elif dimensions == '3D':
        fig = pp.figure(figsize = figsize)
        ax = fig.add_subplot(111, projection='3d', proj_type = proj)
    else:
        print("Invalid dimensions, supports 2D or 3D")
    return fig,ax

def DressFig(fig, ax, title='', savefig = False, folder = '', name = 'placeholder.png',
                 xlabel = None, ylabel = None , zlabel=None, xlim = None, ylim = None,  zlim = None,
                 legendlist = None, logScale = False, logScaleX = False):
    """Dress up a figure and axes set with various commands to adjust its styling"""
    ax.set_title(title)    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if zlabel:
        ax.set_zlabel(zlabel)
    
    if logScale:
        ax.set_yscale('log')
    if logScaleX:
        ax.set_xscale('log')
    
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if zlim:
        ax.set_zlim(zlim)
    
    if legendlist != None:
        ax.legend(legendlist)
    
    if savefig: 
        saveFig(folder,name)
    
    return fig,ax

def Animate3DPlotToGIF(fig, ax, filename, folder):
    """Creates a /GIF/ folder in the given dir to store the frames and results
    
    Depricated method, preferred approach is class based now"""
    
    folder = folder + 'GIF/'
    folder_frames = folder + 'Frames/'
    
    LL.CheckDir(folder)
    LL.CheckDir(folder_frames)
    
    from socket import gethostname
    computerID = gethostname();
    
    if computerID == 'DESKTOP-2477IF6': # Laptop
        ffmpeg_loc = "C:/Users/mrfan/Anaconda3/pkgs/ffmpeg-4.2.2-he774522_0/Library/bin/"
    elif computerID == 'DESKTOP-VDTVGEC': # desktop
        ffmpeg_loc = "D:/Anacond/pkgs/ffmpeg-4.2-ha925a31_0/Library/bin/"
    else:
        print("THis computer is not set up for animating a 3D plot!")
        return
    
    input_loc = folder_frames
    output_loc = folder
    
    for angle in range(0,360,2):
        i = np.int((angle)/2)
        ax.view_init(30,angle)
    
        framename= "frame" + str(i) + ".jpg"
        pp.savefig(folder_frames + framename, dpi=96)
        pp.gca()
    
    
    
    command = ffmpeg_loc + "ffmpeg -i \""+input_loc+"frame%d.jpg\" \""+ output_loc + filename + ".gif\""
    command.replace("//", "/")
    command.replace(r"/", r"\\")
    os.system(command)
    
    print("Command was: " + command)
    print( "Wrote gif file to: "+ output_loc + filename + ".gif")
    
    return

#%% Add widgets
def AddAxesWidget(ax, pos=None, size=1):
    """Add a 3d widget to a 3-d figure showing the directions of the axes"""
    if pos == None:
        pos = [0,0,0]
        
    
    ax.quiver(pos[0], pos[1], pos[2], 1,0,0, length=size)
    ax.quiver(pos[0], pos[1], pos[2], 0,1,0, length=size)
    ax.quiver(pos[0], pos[1], pos[2], 0,0,1, length=size)
    
    return

def AddVector(ax,p1,p2, color=None):
    """Add a single vector between p1 and p2 to axes ax"""
    ax.quiver(p1[0],p1[1],p1[2],p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2], length=1, color=color)
    return


