# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:49:30 2020

@author: USER
"""
# import csv
import numpy as np
import pandas as pd
def ImportPropagationHistory(file, ignoreFinalCols = 0, exportMoon = False):
    # import a general TUDAT propagation history file
    

    with open(file, 'r') as file:                
        
         df= pd.read_csv(file,header = None)
         df = df.rename(columns={x:y for x,y in zip(df.columns,range(0,len(df.columns)))})
         t = df.loc[:,0].values
         x = df.loc[:,1::6].values
         y = df.loc[:,2::6].values
         z = df.loc[:,3::6].values
         vx = df.loc[:,4::6].values
         vy = df.loc[:,5::6].values
         vz = df.loc[:,6::6].values
         moonState = []
                   
         if ignoreFinalCols !=0:
             if exportMoon:
                 moonState = np.array([ x[:,-1], y[:,-1], z[:,-1], vx[:,-1], vy[:,-1], vz[:,-1]  ])
                 moonState = moonState.T
                 
             x = x[:,:-ignoreFinalCols];
             y = y[:,:-ignoreFinalCols];
             z = z[:,:-ignoreFinalCols];
             vx = vx[:,:-ignoreFinalCols];
             vy = vy[:,:-ignoreFinalCols];
             vz = vz[:,:-ignoreFinalCols];
             
             
    return t, x, y, z, vx, vy, vz ,moonState

def ImportMoonHistory(file):
    # import a general TUDAT propagation history file
    

    with open(file, 'r') as file:                
        
         df= pd.read_csv(file,header = None)
         df = df.rename(columns={x:y for x,y in zip(df.columns,range(0,len(df.columns)))})
         x = df.loc[:,1].values
         y = df.loc[:,2].values
         z = df.loc[:,3].values
         vx = df.loc[:,4].values
         vy = df.loc[:,5].values
         vz = df.loc[:,6].values
         moonState = np.array([ x, y, z, vx, vy, vz  ])
         moonState = moonState.T
             
             
    return moonState

def ImportSingleColumnFile(file):
    # import a single column file
    

    with open(file, 'r') as file:                
        
         df= pd.read_csv(file, header = None)
         df = df.rename(columns={x:y for x,y in zip(df.columns,range(0,len(df.columns)))})
         x = df.loc[:,0].values

             
             
    return x

def ImportDataFile(file, delim = '\t'):
    # import a general TUDAT propagation history file

    with open(file, 'r') as file:                
        
         df= pd.read_csv(file, header = None, sep = delim )
         df = df.rename(columns={x:y for x,y in zip(df.columns,range(0,len(df.columns)))})
         x = df.loc[:,1:].values

             
             
    return x
    
def ConvertToRelativeFrame (f_x, f_y, f_z, x, y, z):
    N = x.shape[0]
    W = x.shape[1]
    x = x - f_x.reshape(N,1) *np.ones((N,W))
    y = y - f_y.reshape(N,1) *np.ones((N,W))
    z = z - f_z.reshape(N,1) *np.ones((N,W))
    
    
    return x,y,z
    
def ConvertToBarycentricFrame(body2_p, body3_p, body2_m = 7.34767309E22, body1_p = None, body1_m =5.97237E24, exportL4 = False, fixedL4 = False):
    """ Convert a 3 dimensional set of body coordinates to a barycentric frame between body 1 and body 2 \n
    body2_p: Cartesian coordinates of the secondary body (moon) \n
    body3_p: Cartesian coordinates of the tertiary body (satellites), multiple satellites supported as [x1, x2, y1,y2,z1,z2] \n
    body2_m: Mass of secondary body, default Lunar mass \n
    body1_p: Position vector of primary body, None uses a centralized form (e.g. Earth centralized frame) \n
    body1_m: Mass of primary body, default Earth mass \n
    exportL4: Export the coordinates of the L4 point in the barycentric frame."""
    
    N = np.int(body3_p.shape[1]/3)
    
    x = body3_p[:,0:N]
    y = body3_p[:,N:2*N]
    z = body3_p[:,2*N:3*N]
    
    if np.any(body1_p == None):
        body1_p = np.zeros(body2_p.shape)
        
    xe, ye, ze = (np.zeros((x.shape)) for i in range(0,3))
    
    if len(x.shape) > 1:
        N = x.shape[0]
        W = x.shape[1]
    else:
        N = x.shape[0]
        W = 1
    
    H = np.cross(body2_p[:,:3], body2_p[:,3:], axis = 1)
    H /= np.linalg.norm(H, axis=1).reshape(N,1)
    
    mu = (body2_m)/(body1_m + body2_m) # 0.0112#0.0385209
    barycentre =  mu*(body2_p[:,:3]-body1_p[:,:3]) + body1_p[:,:3]
    

    x_dir = (barycentre)/np.linalg.norm(barycentre,axis=1).reshape(N,1)
    z_dir = H
    y_dir = np.cross(z_dir,x_dir, axis=1)
    y_dir /= np.linalg.norm(y_dir, axis=1).reshape(N,1)
    
    rel_x = x - barycentre[:,0].reshape(N,1)*np.ones((N,W))
    rel_y = y - barycentre[:,1].reshape(N,1)*np.ones((N,W))
    rel_z = z - barycentre[:,2].reshape(N,1)*np.ones((N,W))
    
    rel = np.concatenate([rel_x,rel_y,rel_z], axis=1)
    
    if W > 1:
        for i in range(0,W):
            state = np.concatenate([rel[:,0+i].reshape(N,1), rel[:,W+i].reshape(N,1),
                                    rel[:,2*W+i].reshape(N,1)], axis=1)       
            for j in range(0,N):
                xe[j,i] = np.dot(state[j,:],x_dir[j,:].T)
                ye[j,i] = np.dot(state[j,:],y_dir[j,:].T)
                ze[j,i] = np.dot(state[j,:],z_dir[j,:].T)
    else:
        state = np.concatenate([rel[:,0].reshape(N,1), rel[:,1].reshape(N,1),
                                    rel[:,2].reshape(N,1)], axis=1)       
        for j in range(0,N):
            xe[j] = np.dot(state[j,:],x_dir[j,:].T)
            ye[j] = np.dot(state[j,:],y_dir[j,:].T)
            ze[j] = np.dot(state[j,:],z_dir[j,:].T)
        
        xe = xe.reshape(N,1)
        ye = ye.reshape(N,1)
        ze = ze.reshape(N,1)
    
    L4 = np.array([])
    if exportL4:
        if fixedL4:
            L4 = np.linalg.norm(body2_p[0,:3]*np.ones(body2_p[:,:3].shape),
                                axis=1).reshape(N,1) *  np.array([0.5-mu, 0.5*np.sqrt(3), 0])
        else:
            L4 = np.linalg.norm(body2_p[:,:3], axis=1).reshape(N,1) *  np.array([0.5-mu, 0.5*np.sqrt(3), 0])
        
    return xe,ye,ze, L4
    
def convertFullStateToBarycentric(body2, body3, body2_m = 7.34767309E22, body1 = None, body1_m =5.97237E24, exportL4 = False, fixedL4 = False):
    """ !! NOT VERIFIED YET!! \n
    Convert a set of full state vectors to a barycentric frame between body 1 and body 2 \n
    body2_p: Cartesian coordinates of the secondary body (moon) \n
    body3_p: Cartesian coordinates of the tertiary body (satellites), multiple satellites supported as [x1, x2, y1,y2,z1,z2] \n
    body2_m: Mass of secondary body, default Lunar mass \n
    body1_p: Position vector of primary body, None uses a centralized form (e.g. Earth centralized frame) \n
    body1_m: Mass of primary body, default Earth mass \n
    exportL4: Export the coordinates of the L4 point in the barycentric frame."""
    
    mu = (body2_m)/(body1_m + body2_m)
    
    if np.any(body1 == None):
        body1 = np.zeros(body2.shape)
        
    b1_pos = body1[:,:3]
    b1_vel = body1[:,3:6]
    
    W = body3.shape[1]
    N = body3.shape[0]
    
    b2_pos = body2[:,:3]
    b2_vel = body2[:,3:6]
    
    b3_pos = body3[:,:3]
    b3_vel = body3[:,3:6] if (W==6) else []

    H = np.cross(b2_pos, b2_vel, axis = 1)
    H /= np.linalg.norm(H, axis=1).reshape(N,1)  
    
    
    # position of the barycentre
    bc =  mu*(b2_pos-b1_pos) + b1_pos
    
    # Grab angular speed
    r = np.linalg.norm((b2_pos),axis=1) **2
    r = r.reshape(N,1)
    r = np.concatenate((r,r,r),axis=1)
    omega = np.cross(b2_pos ,b2_vel ) / r
    
    # velocity of the barycentre
    bcv = np.cross(omega,bc)

    x_dir = (bc)/np.linalg.norm(bc,axis=1).reshape(N,1)
    z_dir = H
    y_dir = np.cross(z_dir,x_dir, axis=1)
    y_dir /= np.linalg.norm(y_dir, axis=1).reshape(N,1)
    
    # Compute relative velocities
    rel_pos = b3_pos - bc
    x, y, z = (np.zeros((N,1)) for i in range(0,3))
    for j in range(0,N):
            x[j] = np.dot(rel_pos[j,:],x_dir[j,:].T)
            y[j] = np.dot(rel_pos[j,:],y_dir[j,:].T)
            z[j] = np.dot(rel_pos[j,:],z_dir[j,:].T)
    pos = np.concatenate((x,y,z), axis = 1)
    
    
    if b3_vel != []:
        rel_vel = b3_vel - np.cross(omega,b3_pos)# compensate for rotating ref. frame
        for j in range(0,N):
            x[j] = np.dot(rel_vel[j,:],x_dir[j,:].T)
            y[j] = np.dot(rel_vel[j,:],y_dir[j,:].T)
            z[j] = np.dot(rel_vel[j,:],z_dir[j,:].T)
            # compensate for rotation of frame itself
        vel = np.concatenate((x,y,z), axis = 1)
    else :
        vel = []
    
    return pos, vel

def computeL4Location(body2, body1 = None, mu = 0.0123):
    N = body2.shape[0]
    if body1 == None: 
        body1 = np.zeros(body2.shape)
    
    p = body2 - body1
    
    M = np.cross(body2[:,:3], body2[:,3:6], axis=1)
    x_dir = (p[:,:3])/np.linalg.norm(p[:,:3],axis=1).reshape(N,1)
    z_dir = M/np.linalg.norm(M)
    y_dir = np.cross(z_dir,x_dir, axis=1)
    y_dir /= np.linalg.norm(y_dir, axis=1).reshape(N,1)
    sc = (np.linalg.norm(p[:,:3], axis=1).reshape(N,1))
    L4 = body1[:,:3] + sc*x_dir*(0.5)*np.ones(p[:,:3].shape) + sc*y_dir*0.5*np.sqrt(3)*np.ones(p[:,:3].shape) 
    return L4