# -*- coding: utf-8 -*-
"""
General call function for the verification of other algorithms

Will attempt to verify the correct working of the defined subfunctions, and 
report its results

Created on Wed Mar 11 10:32:27 2020

@author: Jurriaan van 't Hoff
"""


import VVFunctions as VV


failureCount = 0;
failures = []

#%% DATA IMPORT FUNCTIONS
#%% Relative frame conversion
name = "Relative frame conversion"
verified = VV.VerifyRelativeFrame()

if verified:
    print("%s succesfully validated!" %name)
else:
    print("%s failed validation!" %name)
    failures.append(name)

#%% Baricentric frame conversion
name = "Barycentric frame conversion"
verified = VV.VerifyBarycentric()

if verified:
    print("%s succesfully validated!" %name)
else:
    print("%s failed validation!" %name)
    failures.append(name)

    
#%% GENERAL FUNCTIONS
#%% Get Baselines
name = "Baseline retrieval"
verified = VV.VerifyBaselines()

if verified:
    print("%s succesfully validated!" %name)
else:
    print("%s failed validation!" %name)
    failures.append(name)

#%% 2 Dimensional plane projection
name = "2D plane projection"
verified = VV.Verify2DProj()

if verified:
    print("%s succesfully validated!" %name)
else:
    print("%s failed validation!" %name)
    failures.append(name)

#%% PSF Analysis
name = "PSF Computation"
verified = VV.VerifyPSF()

if verified:
    print("%s succesfully validated!" %name)
else:
    print("%s failed validation!" %name)
    failures.append(name)

#%% PSF Analysis - Slow
#   DEPRECATED FUNCTION

#name = "PSF Computation - Slow"
#verified = VV.VerifyPSFSlow()
#
#if verified:
#    print("%s succesfully validated!" %name)
#else:
#    print("%s failed validation!" %name)
#    failures.append(name)
#    
    
#%% L4 Coordinate conversion
name = "L4 coordinate computation"
verified = VV.VerifyL4Coordinates()

if verified:
    print("%s succesfully validated!" %name)
else:
    print("%s failed validation!" %name)
    failures.append(name)


    
#%% Tally up the score    
    
failureCount = len(failures)

if failureCount == 0:
    print("All functions were succesfull validated, succes!")
else:
    print("Validation found %i total function errors" %failureCount)
    print("The following functions failed verification:")
    print(failures)
    