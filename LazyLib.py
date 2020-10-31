# -*- coding: utf-8 -*-
"""
Lazy Library - A library of quality of life functions that automate annoying tasks.

Created on Thu Feb 27 11:17:24 2020

@author: Jurriaan van 't Hoff
"""

import os, socket
from shutil import copy

def SyncDataFiles(targetedTransfer = False, origDir = None, targetDir = None, copyTargetDirContents = False):
    """    Synchronize TUDAT data results with the local Python data folder \n
    Default operation only synchronizes tudat data output files to the python local folders \n
    Can be used as a general folder copy+paste automisation using targetedTransfer"""
    
    if(targetedTransfer == False):
        computerID = socket.gethostname();
        
        if targetDir == None:
            targetDir = os.getcwd() + '/Data/'

        if computerID == 'DESKTOP-2477IF6': # Laptop
            origDir = 'C:/tudat/tudatBundle/tudatApplications/TestApp/SimulationOutput/'
        elif computerID == 'DESKTOP-VDTVGEC': # desktop
            origDir = 'D:/TUDAT/tudatBundle/tudatApplications/TestApp/SimulationOutput/'
        else:
            return
            
                
    contents = os.listdir(origDir)
    subdirs = [i for i in contents if '.' not in i]
    
    if(copyTargetDirContents):
        CheckDir(targetDir)
        for file in contents: 
            copy(origDir + '/' + file, targetDir)
    
    for directory in subdirs:
        newDataDir = targetDir + directory
        CheckDir(newDataDir)
        filestocopy = os.listdir(origDir + directory)
        filespresent = os.listdir(newDataDir)
        
        for file in filestocopy:
            if ("nosync.txt" in filespresent):
                break;
                
            if (file not in filespresent or 
                os.path.getsize(origDir + directory+ '/'  + file) != 
                os.path.getsize(newDataDir + '/' + file)):
                
                copy(origDir + directory + '/' + file , newDataDir)       
    return
    
def SyncDataOptimization(targetedTransfer = False, origDir = None, targetDir = None, copyTargetDirContents = False):
    """    Synchronize TUDAT data results with the local Python data folder \n
    Default operation only synchronizes tudat data output files to the python local folders \n
    Can be used as a general folder copy+paste automisation using targetedTransfer"""
    
    if(targetedTransfer == False):
        computerID = socket.gethostname();
        
        if targetDir == None:
            targetDir = os.getcwd() + '/Data/'

        if computerID == 'DESKTOP-2477IF6': # Laptop
            origDir = 'C:/tudat/tudatBundle/tudatApplications/Optimization/SimulationOutput/'
        elif computerID == 'DESKTOP-VDTVGEC': # desktop
            origDir = 'D:/TUDAT/tudatBundle/tudatApplications/Optimization/SimulationOutput/'
        else:
            return
            
                
    contents = os.listdir(origDir)
    subdirs = [i for i in contents if '.' not in i]
    
    if(copyTargetDirContents):
        CheckDir(targetDir)
        for file in contents: 
            copy(origDir + '/' + file, targetDir)
    
    for directory in subdirs:
        newDataDir = targetDir + directory
        CheckDir(newDataDir)
        filestocopy = os.listdir(origDir + directory)
        filespresent = os.listdir(newDataDir)
        
        for file in filestocopy:
            if ("nosync.txt" in filespresent):
                break;
                
            if (file not in filespresent or 
                os.path.getsize(origDir + directory+ '/'  + file) != 
                os.path.getsize(newDataDir + '/' + file)):
                
                copy(origDir + directory + '/' + file , newDataDir)       
    return

def CheckDir(d):
    """Check if a directory exists, if not create a new one"""
    
    if not os.path.exists(d):
        os.makedirs(d)
        print(f'Directory did not exist, new one created: {d}')
    return

def RelativeFolder(d):
    """Yields a string to a folder pathed relative to the cwd, and creates it if necessary"""
    
    for i in range(0,d.count('../')):
        os.chdir("../")
        d = d[3:]
        # e = d[3*i:]
        
    cwd = os.getcwd()
    relfol = cwd + d
    CheckDir(relfol)
    
    return relfol

def UpdateFiguresToFolders():
    
    
    return

def CreateFigureDirs():
    for i in range(0,len(Names)):
        D = folder_dict[Names[i]]
        CheckDir(D)    
    return


owd = os.getcwd()

if owd[-6:] != 'Python':
    os.chdir('../../')
else:
    os.chdir('../')

os.chdir('PRESENTATION/media')

pres_media_dir = os.getcwd()

os.chdir('../../Python')

if owd[-6:] != 'Python':
    os.chdir('../')

os.chdir('../Thesis paper/Figures/')
    
d = os.getcwd()
folder_dict ={
    "Abstract": d+ "/Abstract/",
    "Title": d+"/0_Title/",
    "Preface": d + "/1_Preface/",
    "Introduction": d +"/2_Introduction/",
    "RadioAstronomy": d + "/3_RadioAstronomy/",
    "OLFAR": d + "/4_OLFAR/",
    "SARI": d + "/5_SARI/",
    "Science": d + "/6_Science/",
    "MissionRequirements": d + "/7_MissionRequirements/",
    "Deployment": d + "/8_Deployment/",
    "Methodology": d + "/9_Methodology/",
    "NumericalSim": d + "/10_NumericalSim/",
    "Results": d + "/11_Results/",
    "Conclusion": d + "/12_Conclusion/",
    "Appendix": d + "/A_Appendix/",
    "Pres": d + "/Presentation/",
    "Pres_media": pres_media_dir}

d = owd + "/Figures/"
folder_dict_funcs ={
    "PA": d + "/Perturbation Analysis/",
    "SD": d + "/SampleDistribution/",
    "VER": d + "/Verification/",
    "BP": d + "/BodyPointing/",
    "IA": d + "/Integrator Analysis/",
    "RK4": d + "/RK4 Analysis/",
    "CF": d + "/Cost Function/",
    "CM": d + "/Cost Mapping/",
    "MISC": d + "/MISC/",
    "OPT": d + "/Optimization Results/"}


Names = ["Preface", "Title", "Preface", "Introduction", "RadioAstronomy",
        "OLFAR","SARI","Science","MissionRequirements", "Deployment",
        "Methodology","NumericalSim","Results","Conclusion",
        "Appendix"]

Constants ={
    "c": 	299792458,
    "Lunar_mean_dist": 385000e3,
    "Lunar_semi_major_axis": 384748e3,
    "Lunar_mass": 7.342e22,
    "mu*": 0.012300,
    "mu": 0.9877,
    "AU": 1.495978707E11,
    "G": 6.6743015E-11
    }
#CreateFigureDirs()

os.chdir(owd)

# SyncDataFiles()
