# -*- coding: utf-8 -*-
"""
Container file with code for class-based figures for easy animation generation

Figures can be created as a derived subclass from modularfigure with function drawuntilentry()
defining the general plotting script for the figure. These figures can then be animated using the defined
AnimateXXX functions

Created on Mon Oct  5 14:15:02 2020

@author: Jurriaan
"""

import numpy as np
import sys 

import os
from os import chdi

import subprocess

from matplotlib import pyplot as pp
from mpl_toolkits.mplot3d import Axes3D # Needs to be imported to support older python versions
from matplotlib.ticker import FormatStrFormatter

from tqdm import tqdm

import pandas as pd

sys.path.append('../')
chdir('../')
from LazyLib import CheckDir
from LazyLib import Constants
import PlottingFunctions as PF

chdir('Presentation Prep')

#%% Support functions


def vectorProj(v1,v2):
    """Projection of vector v1 onto v2"""
    
    v3 = np.dot(v1,v2)/np.linalg.norm(v2)
    return v3



#%% Animate plot

def GrabFolders(folder):
    """Create folder for the output of the animation and its frames, and retrieve
    the location of ffmpeg.exe based on the socket of the computer"""
    folder_frames = folder + 'Frames/'
    
    CheckDir(folder)
    CheckDir(folder_frames)
    
    # Find the install of FFPMEG required to render gifs
    from socket import gethostname
    computerID = gethostname();
    
    ########
    # Code necessary to assist with working from different locations
    # Computer's unique socket ID is tied to the location of ffmpeg
    ########
    if computerID == 'DESKTOP-2477IF6': # Laptop
        ffmpeg_loc = "C:/Users/mrfan/Anaconda3/pkgs/ffmpeg-4.2.2-he774522_0/Library/bin/"
    elif computerID == 'DESKTOP-VDTVGEC': # desktop
        ffmpeg_loc = "D:/ProgramData/Anaconda3/pkgs/ffmpeg-4.2.2-he774522_0/Library/bin/"
    else:
        print("This computer is not set up for animating a 3D plot!")
        return
    
    

    return folder_frames, folder, ffmpeg_loc

def AnimateRotating3DPlot(fig, ax, filename, folder):
    """Creates a /GIF/ folder in the given dir to store the frames and results"""
    
    input_loc, output_loc, ffmpeg_loc = GrabFolders(folder)
    
    print("Generating .gif frames: \n")
    
    for angle in tqdm(range(0,360,2)):
        i = np.int((angle)/2)
        ax.view_init(30,angle)
    
        framename= "frame" + str(i) + ".jpg"
        pp.savefig(input_loc + framename, dpi=96)
        pp.gca()
    
    
    
    command = ffmpeg_loc + "ffmpeg -y -i \""+input_loc+"frame%d.jpg\" \""+ output_loc + filename + ".gif\""
    command.replace("//", "/")
    command.replace(r"/", r"\\")
    os.system(command)
    
    print("Command was: " + command)
    print( "Wrote gif file to: "+ output_loc + filename + ".gif")
    
    return

def AnimateDataRangePlot(modularfig, filename, folder,rnge = 1e8,skip = 1):
    """Creates a /GIF/ folder in the given dir to store the frames and results"""
    
    input_loc, output_loc, ffmpeg_loc = GrabFolders(folder)
    
    # Interact with mod. figure
    
    # retrieve iteratable range
    

    it = min(rnge,modularfig.iterablerange)
        
        
    modularfig.outputfolder = input_loc
    
    print("Generating .gif frames: \n")
    
    for i in tqdm(range(0,int(it/skip))):
        framename= "frame" + str(i) + ".jpg"
        
        modularfig.plotUntilEntry(i*skip)
        modularfig.filename = framename
        
        modularfig.saveFigure()
        modularfig.close()
        
    
    
    
    command = ffmpeg_loc + "ffmpeg -y -i \""+input_loc+"frame%d.jpg\" \""+ output_loc + filename + ".gif\""
    command.replace("//", "/")
    command.replace("/", r"\\")
    
    subprocess.call(command)    
    
    command = ffmpeg_loc + "ffmpeg -y -i \""+input_loc+"frame%d.jpg\" \""+ output_loc + filename + ".mp4\""
    command.replace("//", "/")
    command.replace("/", r"\\")
    
    subprocess.call(command)
    
    # print("Command was: " + command)
    print( "Wrote files to: "+ output_loc + filename + ".gif")
    
    return

#%% Class structure for plots

class modularfigure:
    
    
    
    def __init__(self, fig, filename = 'placeholder', outputfolder = '' ):
        self.figure = fig
        self.filename = filename
        self.outputfolder = outputfolder
        self.iterablerange = 1
        
    def setOutput(self, filename,outputfolder):
        self.filename = filename
        self.outputfolder = outputfolder
    def saveFigure(self):
        self.figure.savefig(self.outputfolder + self.filename)
        
    def show(self):
        pp.show()
    def close(self):
        pp.close(self.figure)
        
        
class combinedOrbitFigure(modularfigure):
    #Data arrays

    
    
    def __init__(self):
        
        
        self.figure, self.axes = pp.subplots(1,3, figsize = (15,5))
        
        super(combinedOrbitFigure,self).__init__(self.figure)
        self.x , self.y, self.z = [],[],[]
        self.xlim, self.ylim, self.zlim = [-.4, 1.2], [-.1,1.3], [-2, 2]
        self.scl = Constants["Lunar_mean_dist"];
        self.mu = Constants["mu*"]
        self.arr = []
        
    def updateData (self, x,y,z, dcm=1):

        self.x, self.y, self.z = x[::dcm],y[::dcm],z[::dcm]
        self.iterablerange = len(self.x)
        
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltx, self.plty, self.pltz = \
            np.nan*np.empty_like(self.x),np.nan*np.empty_like(self.y),np.nan*np.empty_like(self.z)
        self.setupFigure()
        
        self.times = pd.date_range('2030-01-01', periods=len(self.x), freq='240min')
        
    def setupFigure(self):
        """Setup the figure before plotting """
        
        self.lnxy, = self.axes[0].plot(np.nan*self.x,np.nan*self.y)
        self.mkxy = self.axes[0].scatter(0,0,c='r', marker='s',zorder=11)
        self.axes[0].scatter(-self.mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
        self.axes[0].scatter(1-self.mu,0,marker = 'o',c = 'g', zorder=10)
        self.axes[0].scatter(0.5,np.sqrt(3)/2,marker = 'd', c = 'k', zorder=10)
        self.axes[0].set_xlabel('$x_B$ [R]')
        self.axes[0].set_ylabel('$y_B$ [R]')
        self.axes[0].set_xlim(self.xlim)
        self.axes[0].set_ylim(self.ylim)
        
        self.lnyz, = self.axes[2].plot(np.nan*self.x,np.nan*self.y)
        self.mkyz = self.axes[2].scatter(0,0,c='r', marker='s',zorder=11)
        self.axes[2].scatter(0,0,marker = 'o', s = 70, c = 'r', zorder=10)
        self.axes[2].scatter(0,0,marker = 'o',c = 'g', zorder=10)
        self.axes[2].scatter(np.sqrt(3)/2,0,marker = 'd', c = 'k', zorder=10)
        self.axes[2].set_xlabel('$y_B$ [R]')
        self.axes[2].set_ylabel('$z_B$ [R * 1E-2]')
        self.axes[2].set_xlim(self.ylim)
        self.axes[2].set_ylim(self.zlim)
        self.axes[2].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        self.lnxz, = self.axes[1].plot(np.nan*self.x,np.nan*self.y)
        self.mkxz = self.axes[1].scatter(0,0,c='r', marker='s',zorder=11)
        self.axes[1].scatter(-self.mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
        self.axes[1].scatter(1-self.mu,0,marker = 'o',c = 'g', zorder=10)
        self.axes[1].scatter(0.5,0,marker = 'd', c = 'k', zorder=10)
        self.axes[1].set_xlabel('$x_B$ [R]')
        self.axes[1].set_ylabel('$z_B$ [R * 1E-2]')
        self.axes[1].set_xlim(self.xlim)
        self.axes[1].set_ylim(self.zlim)
        self.axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        
        
        self.datetext = self.axes[1].text(0.5,-1.8,'')
        self.figure.legend(['Swarm \n motion', 'Swarm', 'Earth', 'Moon', 'L4' ],loc=(0.72,0.67))
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.axes[0].set_xlim(self.xlim)
        self.axes[0].set_ylim(self.ylim)
        self.axes[1].set_xlim(self.xlim)
        self.axes[1].set_ylim(self.zlim)
        self.axes[2].set_xlim(self.ylim)
        self.axes[2].set_ylim(self.zlim)

    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltx[:i] = self.x[:i]/self.scl
        self.plty[:i] = self.y[:i]/self.scl
        self.pltz[:i] = self.z[:i]*100/self.scl

        #Update data of lines and markers
        self.lnxy.set_xdata(self.pltx) 
        self.lnxy.set_ydata(self.plty) 
         
        self.arr = [[self.x[i][0]/self.scl , self.y[i][0]/self.scl]]
        self.mkxy.set_offsets(self.arr)
       
        self.lnyz.set_xdata(self.plty) 
        self.lnyz.set_ydata(self.pltz) 
        
        self.arr = [[self.y[i][0]/self.scl , self.z[i][0]*100/self.scl]]
        self.mkyz.set_offsets(self.arr)
        
        self.lnxz.set_xdata(self.pltx) 
        self.lnxz.set_ydata(self.pltz) 
        
        self.arr = [[self.x[i][0]/self.scl , self.z[i][0]*100/self.scl]]
        self.mkxz.set_offsets(self.arr)
        
        #Update to the proper datetime
        self.datetext.set_text(self.times[i])


class OrbitFigure(modularfigure):
    #Data arrays

    
    
    def __init__(self):
        
        self.figure1, self.axes0 = pp.subplots(1,1, figsize = (5,5))
        # self.figure2, self.axes1 = pp.subplots(1,1, figsize = (5,5))
        # self.figure3, self.axes2 = pp.subplots(1,1, figsize = (5,5))
        
        super(OrbitFigure,self).__init__(self.figure1)
        self.x , self.y = [],[]
        self.xlim, self.ylim, self.zlim = [-.4, 1.2], [-.1,1.3], [-2, 2]
        self.scl = Constants["Lunar_mean_dist"];
        self.mu = Constants["mu*"]
        self.arr = []
        
        
    def saveFigure(self):
        self.figure1.savefig(self.outputfolder + self.filename)
        # self.figure2.savefig(self.outputfolder + self.filename)
        # self.figure3.savefig(self.outputfolder + self.filename)
        
    def updateData (self, x,y,z, dcm=1):

        self.x, self.y, self.z = x[::dcm],y[::dcm],z[::dcm]
        self.iterablerange = len(self.x)
        
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltx, self.plty, self.pltz = \
            np.nan*np.empty_like(self.x),np.nan*np.empty_like(self.y),np.nan*np.empty_like(self.z)
        self.setupFigure()
        
        self.times = pd.date_range('2030-01-01', periods=len(self.x), freq='240min')
        
    def setupFigure(self):
        """Setup the figure before plotting """
        
        self.lnxy, = self.axes0.plot(np.nan*self.x,np.nan*self.y)
        self.mkxy = self.axes0.scatter(0,0,c='r', marker='s',zorder=11)
        self.axes0.scatter(-self.mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
        self.axes0.scatter(1-self.mu,0,marker = 'o',c = 'g', zorder=10)
        self.axes0.scatter(0.5,np.sqrt(3)/2,marker = 'd', c = 'k', zorder=10)
        self.axes0.set_xlabel('$x_B$ [R]')
        self.axes0.set_ylabel('$y_B$ [R]')
        self.axes0.set_xlim(self.xlim)
        self.axes0.set_ylim(self.ylim)
        
        # self.lnyz, = self.axes2.plot(np.nan*self.x,np.nan*self.y)
        # self.mkyz = self.axes2.scatter(0,0,c='r', marker='s',zorder=11)
        # self.axes2.scatter(0,0,marker = 'o', s = 70, c = 'r', zorder=10)
        # self.axes2.scatter(0,0,marker = 'o',c = 'g', zorder=10)
        # self.axes2.scatter(np.sqrt(3)/2,0,marker = 'd', c = 'k', zorder=10)
        # self.axes2.set_xlabel('$y_B$ [R]')
        # self.axes2.set_ylabel('$z_B$ [R * 1E-2]')
        # self.axes2.set_xlim(self.ylim)
        # self.axes2.set_ylim(self.zlim)
        # self.axes2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        # self.lnxz, = self.axes1.plot(np.nan*self.x,np.nan*self.y)
        # self.mkxz = self.axes1.scatter(0,0,c='r', marker='s',zorder=11)
        # self.axes1.scatter(-self.mu,0,marker = 'o', s = 70, c = 'r', zorder=10)
        # self.axes1.scatter(1-self.mu,0,marker = 'o',c = 'g', zorder=10)
        # self.axes1.scatter(0.5,0,marker = 'd', c = 'k', zorder=10)
        # self.axes1.set_xlabel('$x_B$ [R]')
        # self.axes1.set_ylabel('$z_B$ [R * 1E-2]')
        # self.axes1.set_xlim(self.xlim)
        # self.axes1.set_ylim(self.zlim)
        # self.axes1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        
        
        self.datetext = self.axes0.text(0.5,-1.8,'')
        self.figure.legend(['Swarm \n motion', 'Swarm', 'Earth', 'Moon', 'L4' ],loc=(0.72,0.67))
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.axes0.set_xlim(self.xlim)
        self.axes0.set_ylim(self.ylim)


    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltx[:i] = self.x[:i]/self.scl
        self.plty[:i] = self.y[:i]/self.scl
        self.pltz[:i] = self.z[:i]*100/self.scl

        #Update data of lines and markers
        self.lnxy.set_xdata(self.pltx) 
        self.lnxy.set_ydata(self.plty) 
         
        self.arr = [[self.x[i][0]/self.scl , self.y[i][0]/self.scl]]
        self.mkxy.set_offsets(self.arr)
       
        
        #Update to the proper datetime
        self.datetext.set_text(self.times[i])

class relativeToCoreMotion(modularfigure):
    """Plot class for showing swarm motion relative to core as it develops"""

    
    
    def __init__(self):
        
        self.figure = pp.figure(figsize = (6,6))
        self.ax = self.figure.add_subplot(111, projection='3d')

        super(relativeToCoreMotion,self).__init__(self.figure)
        self.x , self.y, self.z = [],[],[]
        self.xlim, self.ylim, self.zlim = [-100, 100], [-100,100], [-100, 100]
        self.scl = 1000
        self.arr = []
        self.lines = []
        self.pl = 0
    def updateData (self, x,y,z, dcm=1):

        self.x, self.y, self.z = x[::dcm,:],y[::dcm,:],z[::dcm,:]
        self.iterablerange = len(self.x)
        self.nsats = self.x.shape[1]
        
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltx, self.plty, self.pltz = \
            np.nan*np.empty_like(self.x),np.nan*np.empty_like(self.y),np.nan*np.empty_like(self.z)
        
        
        self.times = pd.date_range('2030-01-01', periods=len(self.x), freq='240min')
        self.setupFigure()
        
    def setupFigure(self):
        """Setup the figure before plotting """
        
        #create scatter objects
        self.mk = self.ax.scatter(self.x[0],self.y[0],self.z[0],c='r', marker='s',zorder=11,label= 'Satellites')
        self.lines.append(self.ax.plot(np.nan*self.x[0],np.nan*self.y[0],np.nan*self.z[0],c='b',alpha=0.7, label='Trail')[0])

        self.ax.scatter(0,0,0,c= 'k', marker='d', s=80, label='Core')
    
        for i in range(1,self.nsats):
            self.lines.append(self.ax.plot(np.nan*self.x[i],np.nan*self.y[i],np.nan*self.z[i],c='b',alpha=0.7)[0])
            
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        self.ax.set_xlabel('distance from core [km]')
        self.ax.set_ylabel('distance from core [km]')
        self.ax.set_zlabel('distance from core [km]')
        self.title = self.ax.set_title(self.times[0])
        # self.datetext = self.ax.text(80,-95,-95,'s')
        self.figure.legend()
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        
    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltx[:i+1] = self.x[:i+1]/self.scl
        self.plty[:i+1] = self.y[:i+1]/self.scl
        self.pltz[:i+1] = self.z[:i+1]/self.scl
        
        if self.pl > 0 and i > self.pl:
            self.pltx[:i-self.pl] *= np.nan
            self.plty[:i-self.pl] *= np.nan
            self.pltz[:i-self.pl] *= np.nan

        #Update data of lines and markers
        c=0
        for line in self.lines:
            line.set_xdata(self.pltx[:,c])
            line.set_ydata(self.plty[:,c])
            line.set_3d_properties(self.pltz[:,c])
            c += 1
        
        self.mk._offsets3d = (self.pltx[i,:],self.plty[i,:],self.pltz[i,:])
        
        #Update to the proper datetime
        self.title.set_text(self.times[i])
        
class harmonicDecayPlot(modularfigure):
    
    """Plot class for showing swarm motion relative to core as it develops, in 3 stacked plots for coordinates"""

    
    
    def __init__(self):
        
        self.figure, self.axes = pp.subplots(3,1, figsize = (15,7))

        super(harmonicDecayPlot,self).__init__(self.figure)
        self.x , self.y, self.z = [],[],[]
        self.xlim, self.ylim = [0,365], [-50,50]
        self.td = 0
        self.scl = 1000
        self.lines = []

        
    def updateData (self, x,y,z,td, dcm=1):
        """Update data in the plot with support for decimation through dcm"""
        self.td = td[::dcm]
        self.x, self.y, self.z = x[::dcm,:],y[::dcm,:],z[::dcm,:]
        self.iterablerange = len(self.x)
        self.nsats = self.x.shape[1]
        
        self.tdscat = self.td.reshape((len(td),1))*np.ones(self.x.shape)
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltx, self.plty, self.pltz = \
            np.nan*np.empty_like(self.x),np.nan*np.empty_like(self.y),np.nan*np.empty_like(self.z)
        
        
        self.times = pd.date_range('2030-01-01', periods=len(self.x), freq='240min')
        self.setupFigure()
        
    def setupFigure(self):
        """Setup the figure before plotting """
        
    
        for i in range(0,self.nsats):
            self.lines.append(self.axes[0].plot(self.td,np.nan*self.x[:,0],c='b',alpha=0.4)[0])
            self.lines.append(self.axes[1].plot(self.td,np.nan*self.y[:,0],c='b',alpha=0.4)[0])
            self.lines.append(self.axes[2].plot(self.td,np.nan*self.z[:,0],c='b',alpha=0.4)[0])
            
            
        self.axes[0].set_ylabel('Relative $x_B$[km]')
        self.axes[1].set_ylabel('Relative $y_B$[km]')
        
        self.axes[2].set_xlabel('time in orbit [days]')
        self.axes[2].set_ylabel('Relative $z_B$[km]')
        
        self.axes[0].set_xlim(self.xlim)
        self.axes[1].set_xlim(self.xlim)
        self.axes[2].set_xlim(self.xlim)
        
        self.axes[0].set_ylim(self.ylim)
        self.axes[1].set_ylim(self.ylim)
        self.axes[2].set_ylim(self.ylim)
        
        self.mkx = self.axes[0].scatter(self.tdscat[0],self.x[0], c= 'r', marker = 's' )
        self.mky = self.axes[1].scatter(self.tdscat[0],self.y[0], c= 'r', marker = 's')
        self.mkz = self.axes[2].scatter(self.tdscat[0],self.z[0], c= 'r', marker = 's')
        
        self.title = self.axes[0].set_title(self.times[0])
        self.figure.legend(['Satellite motion'])
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        
    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltx[:i+1] = self.x[:i+1]/self.scl
        self.plty[:i+1] = self.y[:i+1]/self.scl
        self.pltz[:i+1] = self.z[:i+1]/self.scl

        # Update plot data            
        for j in range(0,self.nsats):
            self.lines[j*3].set_ydata(self.pltx[:,j])
            self.lines[j*3+1].set_ydata(self.plty[:,j])
            self.lines[j*3+2].set_ydata(self.pltz[:,j])
        
        self.mkx.set_offsets(np.array([self.tdscat[i], self.pltx[i]]).T)
        self.mky.set_offsets(np.array([self.tdscat[i], self.plty[i]]).T)
        self.mkz.set_offsets(np.array([self.tdscat[i], self.pltz[i]]).T)
        
        #Update to the proper datetime
        self.title.set_text(self.times[i])


class uvwBaselinePlot3D(modularfigure):
    """Plot class for showing swarm motion relative to core as it develops"""

    
    
    def __init__(self):
        
        self.figure = pp.figure(figsize = (6,6))
        self.ax = self.figure.add_subplot(111, projection='3d')

        super(uvwBaselinePlot3D,self).__init__(self.figure)
        self.x , self.y, self.z = [],[],[]
        self.xlim, self.ylim, self.zlim = [-100, 100], [-100,100], [-100, 100]
        self.scl = 1000
        self.lines = []
        self.pl = 0
        
    def updateData (self, x,y,z, dcm=1):

        self.x, self.y, self.z = x[::dcm,:],y[::dcm,:],z[::dcm,:]
        self.iterablerange = len(self.x)
        self.nbaselines = self.x.shape[1]
        
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltx, self.plty, self.pltz = \
            np.nan*np.empty_like(self.x),np.nan*np.empty_like(self.y),np.nan*np.empty_like(self.z)
        
        
        self.times = pd.date_range('2030-01-01', periods=len(self.x), freq='240min')
        self.setupFigure()
        
    def setupFigure(self):
        """Setup the figure before plotting """
        
        #create scatter objects
        # self.mk = self.ax.scatter(self.x[0],self.y[0],self.z[0],c='r', marker='s',zorder=11,label= 'Satellites')
        # self.lines.append(self.ax.plot(np.nan*self.x[0],np.nan*self.y[0],np.nan*self.z[0],c='b',alpha=0.7, label='Trail')[0])
        
        self.lines.append(self.ax.plot(np.nan*self.x[:,0],np.nan*self.y[:,0],np.nan*self.z[:,0],c='b',alpha=0.7)[0])
        self.lines.append(self.ax.plot(np.nan*self.x[:,0],np.nan*self.y[:,0],np.nan*self.z[:,0],c='g',alpha=0.7)[0])
        
        for i in range(1,self.nbaselines):
            self.lines.append(self.ax.plot(np.nan*self.x[i],np.nan*self.y[i],np.nan*self.z[i],c='b',alpha=0.7)[0])
            self.lines.append(self.ax.plot(np.nan*self.x[i],np.nan*self.y[i],np.nan*self.z[i],c='g',alpha=0.7)[0])
            
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        self.ax.set_xlabel('[km]')
        self.ax.set_ylabel('[km]')
        self.ax.set_zlabel('[km]')
        self.title = self.ax.set_title(self.times[0])
        # self.datetext = self.ax.text(80,-95,-95,'s')
        self.figure.legend(['Positive baseline',
                            'Negative pair'])
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        
    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltx[:i+1] = self.x[:i+1]/self.scl
        self.plty[:i+1] = self.y[:i+1]/self.scl
        self.pltz[:i+1] = self.z[:i+1]/self.scl
        
        if self.pl > 0 and i > self.pl:
            self.pltx[:i-self.pl] *= np.nan
            self.plty[:i-self.pl] *= np.nan
            self.pltz[:i-self.pl] *= np.nan

        #Update data of lines and markers
        for j in range(0,self.nbaselines):
            # first line is in axes[0], xy plane
            self.lines[j*2].set_xdata(self.pltx[:,j])
            self.lines[j*2].set_ydata(self.plty[:,j])
            self.lines[j*2].set_3d_properties(self.pltz[:,j])
            
            self.lines[j*2+1].set_xdata(-self.pltx[:,j])
            self.lines[j*2+1].set_ydata(-self.plty[:,j])
            self.lines[j*2+1].set_3d_properties(-self.pltz[:,j])
        
        # self.mk._offsets3d = (self.pltx[i,:],self.plty[i,:],self.pltz[i,:])
        
        #Update to the proper datetime
        self.title.set_text(self.times[i])
    
    
    
    
    
class uvwBaselinePlot3Panel(modularfigure):
    def __init__(self):
        
        self.figure, self.axes = pp.subplots(1,3, figsize = (15,5))

        super(uvwBaselinePlot3Panel,self).__init__(self.figure)
        self.x , self.y, self.z = [],[],[]
        self.xlim, self.ylim, self.zlim = [-110, 110], [-110,110], [-110, 110]
        self.scl = 1000
        self.lines = []

        
    def updateData (self, x,y,z, dcm=1):
        """Update data in the plot with support for decimation through dcm"""
        self.x, self.y, self.z = x[::dcm,:],y[::dcm,:],z[::dcm,:]
        self.iterablerange = len(self.x)
        self.nbaselines = self.x.shape[1]
        
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltx, self.plty, self.pltz = \
            np.nan*np.empty_like(self.x),np.nan*np.empty_like(self.y),np.nan*np.empty_like(self.z)
        
        
        self.times = pd.date_range('2030-01-01', periods=len(self.x), freq='240min')
        self.setupFigure()
        
    def setupFigure(self):
        """Setup the figure before plotting """
    
        self.mkxy = self.axes[0].scatter(self.x[0], self.y[0],c='r', marker = 's')
        self.mkyz = self.axes[1].scatter(self.x[0], self.y[0],c='r', marker = 's')
        self.mkxz = self.axes[2].scatter(self.x[0], self.y[0],c='r', marker = 's')
        
        for i in range(0,self.nbaselines):
            self.lines.append(self.axes[0].plot(np.nan*self.x[:,0],np.nan*self.x[:,0],c='b',alpha=0.1)[0])
            self.lines.append(self.axes[1].plot(np.nan*self.x[:,0],np.nan*self.y[:,0],c='b',alpha=0.1)[0])
            self.lines.append(self.axes[2].plot(np.nan*self.x[:,0],np.nan*self.z[:,0],c='b',alpha=0.1)[0])
            
        
        self.axes[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.axes[0].set_xlabel('x [km]')
        self.axes[0].set_ylabel('y [km]')
        self.axes[0].set_xlim(self.xlim)
        self.axes[0].set_ylim(self.ylim)
        
        self.axes[2].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.axes[2].set_xlabel('y [km]')
        self.axes[2].set_ylabel('z [km]')
        self.axes[2].set_xlim(self.ylim)
        self.axes[2].set_ylim(self.zlim)
        
        self.axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.axes[1].set_xlabel('x [km]')
        self.axes[1].set_ylabel('z [km]')
        self.axes[1].set_xlim(self.xlim)
        self.axes[1].set_ylim(self.zlim)
        
        
        self.title = self.axes[1].set_title(self.times[0])
        self.figure.legend(['uvw baselines'])
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        
    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltx[:i+1] = self.x[:i+1]/self.scl
        self.plty[:i+1] = self.y[:i+1]/self.scl
        self.pltz[:i+1] = self.z[:i+1]/self.scl

        # Update plot data            
        for j in range(0,self.nbaselines):
            # first line is in axes[0], xy plane
            self.lines[j*3].set_xdata(self.pltx[:,j])
            self.lines[j*3].set_ydata(self.plty[:,j])
            # second line is axes[1], xz plane
            self.lines[j*3+1].set_xdata(self.pltx[:,j])
            self.lines[j*3+1].set_ydata(self.pltz[:,j])
            # third line is axes[2]. yz plane
            self.lines[j*3+2].set_xdata(self.plty[:,j])
            self.lines[j*3+2].set_ydata(self.pltz[:,j])
            
        self.mkxy.set_offsets(np.array([self.pltx[i], self.plty[i]]).T)
        self.mkyz.set_offsets(np.array([self.plty[i], self.pltz[i]]).T)
        self.mkxz.set_offsets(np.array([self.pltx[i], self.pltz[i]]).T)

        
        #Update to the proper datetime
        self.title.set_text(self.times[i])
        
        
        
class baselineProfile(modularfigure):
    
    def __init__(self):
        
        self.figure, self.ax = PF.MakeFig(figsize=(12,4))
        super(baselineProfile,self).__init__(self.figure)
        self.x , self.y = [],[]
        self.xlim, self.ylim = [0, 1], [0.001,1e9]
        self.scl = 1
        self.lines = []
        self.suplines = []

        
    def updateData (self, R,dR,td, dcm=1):
        """Update data in the plot with support for decimation through dcm"""
        self.td = td[::dcm]
        self.R, self.dR = R[::dcm,:],dR[::dcm,:]
        self.iterablerange = len(self.R[:,0])
        self.nbaselines = int(self.R.shape[1]/2)
        
        self.tdscat = self.td.reshape((len(td),1))*np.ones(self.R.shape)
        
        # create NAN copies of the arrays for plotting
        # to keep all data in the same memory buffer
        self.pltr, self.pltv = \
            np.nan*np.empty_like(self.R),np.nan*np.empty_like(self.dR)        
        
        self.xlim[1] = td[-1]
        self.times = pd.date_range('2030-01-01', periods=len(self.td), freq='240min')
        self.setupFigure()
        
    def setupFigure(self):
        """Setup the figure before plotting """
    
        self.mkr = self.ax.scatter(self.tdscat[0], self.R[0],c='r', marker = 's')
        self.mkv = self.ax.scatter(self.tdscat[0], self.dR[0],c='r', marker = 's')

        
        for i in range(0,self.nbaselines):
            self.lines.append(self.ax.plot(self.td, np.nan*self.R[:,0], color = 'b', alpha = 7.6/self.nbaselines)[0])
            self.lines.append(self.ax.plot(self.td, np.nan*self.dR[:,0], color = 'g', alpha = 7.6/self.nbaselines)[0])
            # self.lines.append(self.axes[0].plot(self.td,np.nan*self.x[:,0],c='b',alpha=0.4)[0])
            # print(f"Added line {i} of {self.nbaselines}")
        # self.lines.append(self.ax.plot(self.td, np.nan*self.R, color = 'b', alpha = 7.6/self.nbaselines))
        # self.lines.append(self.ax.plot(self.td, np.nan*self.dR, color = 'g', alpha = 7.6/self.nbaselines))
        
        self.ax.axhline(500, color = 'r', ls='--', label='Near collision limit')
        self.ax.axhline(100e3, color = 'r',label='100 km drift limit')
        
        self.suplines.append(self.ax.plot(self.td, self.pltr, color = 'k', label='Maximum baseline')[0])
        self.suplines.append(self.ax.plot(self.td, self.pltr, color = 'k', label='Minimum baseline')[0])
        self.suplines.append(self.ax.plot(self.td, self.pltr, color = 'y', label='Mean baseline')[0])
        
        self.suplines.append(self.ax.plot(self.td, self.pltv, color = 'm', label='Maximum baseline rate')[0])
        self.suplines.append(self.ax.plot(self.td, self.pltv, color = 'm', label='Minimum baseline rate')[0])
        self.suplines.append(self.ax.plot(self.td, self.pltv, color = 'm',ls='--', label='Mean baseline rate')[0])
        
        self.ax.axhline(1, color = 'r', label='1 m/s baseline rate')
        self.ax.set_yscale('log')
        
        self.ax.set_xlabel('time since 1 Jan 2030 [days]')
        self.ax.set_ylabel("Baseline magnitude [m] / rate [m/s]")
        
        self.figure.legend( loc=2, bbox_to_anchor=(1.05, 1))
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.title = self.ax.set_title(self.times[0])
        self.figure.tight_layout()   
        
    def updateLims(self, xl, yl, zl):
        self.xlim, self.ylim, self.zlim = xl,yl,zl
        
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_zlim(self.zlim)
        
    def plotUntilEntry(self,i):
        """Render plot to display its data up to a certain integerer entry i in its data arrays. \n
        i must be a integer value within the range of the internal data arrays"""
        # grab more data from the original dataset
        self.pltr[:i+1] = self.R[:i+1]/self.scl
        self.pltv[:i+1] = self.dR[:i+1]/self.scl

        # Update plot data            
        for j in range(0,self.nbaselines):
            # first line is in axes[0], xy plane
            self.lines[j*2].set_ydata(self.pltr[:,j])
            # velocity lines
            self.lines[j*2+1].set_ydata(self.pltv[:,j])

        self.suplines[0].set_ydata(np.max(self.pltr,axis=1))
        self.suplines[1].set_ydata(np.min(self.pltr,axis=1))
        self.suplines[2].set_ydata(np.mean(self.pltr,axis=1))
        
        self.suplines[3].set_ydata(np.max(self.pltv,axis=1))
        self.suplines[4].set_ydata(np.min(self.pltv,axis=1))
        self.suplines[5].set_ydata(np.mean(self.pltv,axis=1))


        self.mkr.set_offsets(np.array([self.tdscat[i], self.pltr[i]]).T)
        self.mkv.set_offsets(np.array([self.tdscat[i], self.pltv[i]]).T)


        
        #Update to the proper datetime
        self.title.set_text(self.times[i])