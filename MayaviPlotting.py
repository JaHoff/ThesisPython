# -*- coding: utf-8 -*-
"""
Plotting function container for all Mayavi related plots

Contains several useful plotting algorithms for different data types
Created on Thu Feb 20 16:17:04 2020

@author: Jurriaan van 't Hoff
"""
import os
import numpy as np
import LazyLib as LL

from mayavi import mlab
from tvtk.api import tvtk # python wrappers for the C++ vtk ecosystem
from tvtk.common import configure_input_data

#%% Figure creation

def MakeFig(figure,size=(600, 600), bgcolor = (0,0,0), fgcolor = (1,1,1), theme = None):
    
    if theme != None:
        colors = themes[theme]
    
        bgcolor = colors[3]
        fgcolor = colors[4]
        
    fig = mlab.figure(figure, bgcolor, fgcolor, size = size)

    return fig

def Title(figure,title):
    mlab.title(title, figure = figure)
    return 
def CloseAll():
    mlab.close(all = True)
    return

#%% Generalized plotting functions

def PlotTextureSphere(fig,texture, R, Rscale=0, x=0, y=0, z=0, Nrad = 180):
    # load and map the texture
    img = tvtk.JPEGReader()
    img.file_name = texture
    texture = tvtk.Texture(input_connection=img.output_port, interpolate=1)
    
    if Rscale !=0:
        R /= Rscale
        x /= Rscale
        y /= Rscale
        z /= Rscale

    # create the sphere source with a given radius and angular resolution
    sphere = tvtk.TexturedSphereSource(radius=R, theta_resolution=Nrad,
                                       phi_resolution=Nrad)

    # assemble rest of the pipeline, assign texture    
    sphere_mapper = tvtk.PolyDataMapper(input_connection=sphere.output_port)
    sphere_actor = tvtk.Actor(mapper=sphere_mapper, texture=texture)
    fig.scene.add_actor(sphere_actor)
    return fig


def PlotColorSphere(fig, Rscale, x=0,y=0,z=0,R = 1, c = (1,1,1), Nrad = 180, opacity = 0.3):
    

    # create the sphere source with a given radius and angular resolution
    sphere = tvtk.SphereSource(center=(x/Rscale, y/Rscale, z/Rscale),radius=R, theta_resolution=Nrad,
                                       phi_resolution=Nrad)

    sphere_mapper2 = tvtk.PolyDataMapper()
    configure_input_data(sphere_mapper2, sphere.output)
    sphere.update()
    p = tvtk.Property(opacity= opacity, color = c)

    actor = tvtk.Actor(mapper=sphere_mapper2, property=p)
    fig.scene.add_actor( actor)
    return sphere

def PlotSurface(fig, x,y,z, lw = 1.0, cmap = 'jet', warp_scale = 1):

    mlab.surf(x,y,z, colormap = cmap, figure = fig, line_width = lw, warp_scale = warp_scale)
    return

def PlotMesh(fig, x,y,z, cmap = 'jet', lw = 2.0, opacity = 1.):
    mlab.mesh(x,y,z,figure = fig, colormap =  cmap, line_width = lw, opacity = opacity)
    return

def Plot3DContour(fig, val, cmap = 'jet', lw = 2.0, opacity = 1., x=[], y=[], z=[]):
    if (len(x) != 0):
        mlab.contour3d(x,y,z,val,figure = fig, colormap =  cmap, line_width = lw, opacity = opacity)
    else:
        mlab.contour3d(val,figure = fig, colormap =  cmap, line_width = lw, opacity = opacity)
    return

def Plot3DQuiver(fig,u,v,w, x=[], y=[], z=[], lw = 3.0, sf = 1, cmap = "jet"):
    if len(x) != 0:
        mlab.quiver3d(x,y,z,u,v,w, figure=fig, line_width = lw, scale_factor = sf,colormap = cmap)
    else:
        mlab.quiver3d(u,v,w, figure=fig, line_width = lw, scale_factor = sf,colormap = cmap)
        
    return

def Plot3D(fig, x, y, z, lw = 3.0, col = (1,1,1), scale = 1000, opacity = 1):
    if len(x.shape) == 2:
        N = x.shape[1]
    else:
        N = 1
    
    plots = dict()
    
    for i in range(0,N):
        plots[i] = mlab.plot3d(x[:,i]/scale, y[:,i]/scale, z[:,i]/scale,
                               color = col,tube_radius=lw, opacity = opacity)
        
    return plots
    
def PlotAnimatedParticles(fig, x,y,z, lw = 3.0, theme = 'default', R =6378e3, speed=1 ):
    """ Plot a satellite swarm with animated orbit progression for a dataset of [H,N] \n
    data points, where N is the number of particles """
    EarthTex = 'EarthTex.jpg'
    
    # switch colors based on theme
    colors = themes[theme]
    
    # Properly scale the data
    xp = x/R 
    yp = y/R
    zp = z/R
    
    L = x.shape[0]
    # Figure out the data format we are dealing with
    if len(x.shape) == 2:
        N = x.shape[1]
    else:
        N = 1
    
    points = dict()
    plots = dict()
    for i in range(0,N):
        points[i] = mlab.points3d(xp[0,i], yp[0,i], zp[0,i],color = colors[2], mode='2ddiamond')
        plots[i] = mlab.plot3d(xp[:,i], yp[:,i], zp[:,i], color = colors[2],tube_radius=lw)
        
        
    # dummy = np.nan * np.empty_like(xp)
    # x2,y2,z2 = np.empty_like(xp), np.empty_like(yp), np.empty_like(zp)
    @mlab.animate(delay=10)
    def anim():
        t = 1
        while True:
            if t >= L: t = 1
            
            # x2[:,:] = dummy[:,:]
            # x2[:t,:] = xp[:t,:]
            # y2[:,:] = dummy[:,:]
            # y2[:t,:] = yp[:t,:]
            # z2[:,:] = dummy[:,:]
            # z2[:t,:] = zp[:t,:]
            
            for j in range(0,N):
                points[j].mlab_source.trait_set(x = xp[t,j], y = yp[t,j], z = zp[t,j])
                # plots[j].mlab_source.trait_set(x = x2[:,j],y = y2[:,j], z = z2[:,j])
                plots[j].mlab_source.reset(x = xp[:t,j], y = yp[:t,j] , z = zp[:t,j])

            t += speed;
            yield
            
    anim()
    return

def AddDefaultAxes(fig = None, color = (1,1,1), Nlabels = 10, lw = 2, xlabel = '', ylabel = '', zlabel = '' ):
    mlab.axes(figure = fig, color = color, nb_labels = Nlabels, line_width = lw,
              xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
    return


#%% Dedicated figure generation functions
def PlotAnimatedEarthOrbit( x,y,z, moon, L4, theme = 'default'):
    """ Plot a satellite swarm with animated orbit progression for a dataset of [H,N] \n
    data points, where N is the number of particles """
    EarthTex = 'EarthTex.jpg'
    
    # switch colors based on theme
    colors = themes[theme]
    
    
    # create a figure window (and scene)
    fig = MakeFig("Satellite swarm orbits around earth",size=(600, 600), 
                  bgcolor = colors[3], fgcolor = colors[4])
    
    
    # Properly scale the data
    R = 6378e3
    xp = x/R 
    yp = y/R
    zp = z/R
    moonp = moon/R
    L4p = L4/R
    
    L = x.shape[0]
    
    # Figure out the data format we are dealing with
    if len(x.shape) == 2:
        N = x.shape[1]
    else:
        N = 1
    
    # Draw Earth
    
    PlotTextureSphere(fig,EarthTex, 1)
    
    # Draw moon and trail (maybe animate)
    moon_trail = mlab.plot3d(moonp[0,0], moonp[0,1], moonp[0,2], color = colors[0])
    moon_dot = PlotColorSphere(fig, 1, moonp[0,0],moonp[0,1],moonp[-1,2], 0.3, c=colors[0] )
    moon_text = mlab.text3d(moonp[0,0], moonp[0,1], moonp[0,2], "Moon", color = colors[0] )
    
   
    # Draw L4 point
    l4 = PlotColorSphere(fig, 1, L4p[-1,0],L4p[-1,1],L4p[-1,2], 0.1, c= colors[1] )
    l4_text = mlab.text3d(L4p[-1,0],L4p[-1,1],L4p[-1,2], "L4", color = colors[1] )
    
    points = dict()
    plots = dict()
    texts = dict()
    for i in range(0,N):
        points[i] = mlab.points3d(xp[0,i], yp[0,i], zp[0,i],color = colors[2],
                                  scale_factor=1, mode = '2ddiamond')
        plots[i] = mlab.plot3d(xp[:1,i], yp[:1,i], zp[:1,i], color = colors[2])
        
    centerpoint = mlab.points3d(xp[0,0], yp[0,0], zp[0,0],color = colors[2],
                                  scale_factor=0.1)
    
    centertext = mlab.text3d(xp[0,0], yp[0,0], zp[0,0], "ARTEMIS", color = colors[2])
        
    @mlab.animate(delay=10)
    def anim():
        t = 1
        while True:
            if t >= L: t = 1
            
            moon_trail.mlab_source.reset( x = moonp[:t,0], y = moonp[:t,1] , 
                             z = moonp[:t,2])
            moon_dot.center = (moonp[t,0],moonp[t,1], moonp[t,2])
            moon_dot.update()
            
            moon_text.position = np.array((moonp[t,0], moonp[t,1], moonp[t,2]))
            # moon_text.mlab_source.x = moonp[t,0]
            # moon_text.mlab_source.trait_set(x=moonp[t,0], y = moonp[t,1], z = moonp[t,2])
            # moon.mlab_source.trait_set(x=moon[t,0], y = moon[t,1], z = moon[t,2])
            
            l4.center = (L4p[t,0], L4p[t,1], L4p[t,2])
            l4.update()
            
            l4_text.position = np.array((L4p[t,0], L4p[t,1], L4p[t,2]))
                
            for j in range(0,N):
                points[j].mlab_source.trait_set(x=xp[t,j], y = yp[t,j], z = zp[t,j])
                plots[j].mlab_source.reset(x = xp[:t,j], y = yp[:t,j] , z = zp[:t,j])
                #texts[j].mlab_source.trait_set(x=xp[t,j], y = yp[t,j], z = zp[t,j])
                
            centerpoint.mlab_source.trait_set(x = xp[t,0], y = yp[t,0] , z = zp[t,0])
            centertext.position = np.array((xp[t,0], yp[t,0] , zp[t,0]))
            t += 1
            yield
            
    anim()


    return fig

def PlotAroundEarth(x,y,z, moondata, L4, plotMoonTrails = False):
    EarthTex = 'EarthTex.jpg'
    # create a figure window (and scene)
    fig = MakeFig("Satellite swarm orbits around earth",size=(600, 600), bgcolor = (0,0,0))
    
    Rearth = 6378e3
    PlotTextureSphere(fig,EarthTex, 1)
    
    width = x.shape[1]
    for i in range(0,width):
        mlab.plot3d(x[:,i]/Rearth,y[:,i]/Rearth,z[:,i]/Rearth)
        PlotColorSphere(fig, Rearth, x[-1,i],y[-1,i],z[-1,i], 0.1, (1,1,0) )
    
    if plotMoonTrails:
        mlab.plot3d(moondata[:,0]/Rearth, moondata[:,1]/Rearth, moondata[:,2]/Rearth, color = (1,0,0))
    
    PlotColorSphere(fig, Rearth, moondata[-1,0],moondata[-1,1],moondata[-1,2], 0.8, c=(1,0,0) )
    PlotColorSphere(fig, Rearth, L4[-1,0],L4[-1,1],L4[-1,2], 0.5, c= (0,0,1) )    


    return fig

#%% Theme setup

# values: mooncolor, l4 color, satcolor, bgcolor, fgcolor
themes = {
    "default": ((1,0,0), (0,0,1), (0,1,0), (0,0,0), (1,1,1)),
    "homeworld": ((1,0,0), (0,1,1), (0,1,0), (0.2,0.2,0.8),(0,1,1))}