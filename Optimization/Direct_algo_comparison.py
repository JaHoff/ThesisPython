# -*- coding: utf-8 -*-
"""
Created on Sun May 31 14:29:34 2020

@author: mrfan
"""


import numpy as np
import sys 

from os import chdir
from os import getcwd

from scipy import spatial


sys.path.append('../')
chdir('../')
import LazyLib as LL
import PlottingFunctions as PF
import DataImport as DI
import GeneralFunctions as GF

figure_folder = LL.folder_dict["NumericalSim"]

#LL.SyncDataOptimization()
print("Synchronized data files")
cwd = getcwd()
datafolder = cwd + '/Data/Optimization/AlgComparison/'

figfolder = cwd + '/Figures/Optimization/'
LL.CheckDir(figfolder)

algorithms = [ "de1220", "sade", "gaco", "pso", "pso_gen"]
algorithms_full = ["Differential Evolution", "Self Adjusting Differential Evolution",
                   "Generational Ant Colony",
                   "Particle Swarm Optimization",
                   "Generational Particle Swarm Optimization"]

n_algos = len(algorithms)
n_gen = 74

j = 0

init = False
for algo in algorithms:
    
    # File = datafolder + 'propagationHistory_'+algo+'_best.dat'
    
    # print("Importing data")
    # t, x, y, z, vx, vy, vz,moon  = DI.ImportPropagationHistory(File,0, False)
    # moon = DI.ImportMoonHistory(datafolder + "propagationHistory_moon.dat")
    # print("Imported data succesfully")
    
    
    
    
    for i in range(1,n_gen):
        file = datafolder + f'fitness_{algo}_sd42_{i}_{i}.dat'
        file_data  = datafolder + f'population_{algo}_sd42_{i}_{i}.dat'
        fitness = DI.ImportSingleColumnFile(file)
        dat = DI.ImportDataFile(file_data)
        
        if init == False: 
            fit = np.zeros((len(fitness),n_gen-1,n_algos))
            fit_mean = np.zeros((1,n_gen-1,n_algos))
            fit_min = np.zeros((1,n_gen-1,n_algos))
            fit_max = np.zeros((1,n_gen-1,n_algos))
            
            pop = np.zeros((dat.shape[0],dat.shape[1], n_algos,n_gen-1))
            resid = np.zeros((n_algos,n_gen-1))
        
        fit[:,i-1,j] = fitness
        
        init=True
    
        pop[:,:,j,i-1] = dat
        
        n_pop = dat.shape[0]
#        simil = np.mean(dat, axis = 0)
#        resid[j,i] = 21/np.sum(np.abs((np.sum(np.abs(dat - simil),axis=0)/(simil*n_pop))))
        
        similarity = 0
        for g in range(1,n_pop):
            similarity += (1 - spatial.distance.cosine(dat[0,:6],dat[g,:6]))
            
        resid[j,i-1] = similarity/n_pop
        
    fit_mean[:,:,j] = np.mean(fit[:,:,j], axis=0)
    fit_max[:,:,j] = np.max(fit[:,:,j],axis=0)
    fit_min[:,:,j] = np.min(fit[:,:,j],axis=0)
    
    
    
    
    
#    fig,ax = PF.MakeFig()
#    PF.Simple2DPlot(gens, fit_mean[:,:,j][0], xmod = 1., ymod = 1., fig = fig, ax = ax)
#    PF.Simple2DPlot(gens, fit_min[:,:,j][0], xmod = 1., ymod = 1., fig = fig, ax = ax)
#    PF.Simple2DPlot(gens, fit_max[:,:,j][0], xmod = 1., ymod = 1., fig = fig, ax = ax)
#    PF.DressFig(fig, ax, title='Optimization progress over generations:' + algo,
#                savefig = True, folder = figfolder, name = 'Direct algo comp generational progress_'+ algo + '.png',
#                     xlabel = "Generation", ylabel = "Cost" , ylim = [1,np.max(fit_max)],
#                     legendlist = ["Mean cost", "Minimal cost", "Maximum cost"],
#                     logScale = False, logScaleX = False)
    
    j += 1

#%%
from matplotlib import pyplot as pp

gens = np.arange(1,n_gen)
cols = ['b', 'r', 'g', 'k', 'y']
fig,ax = PF.MakeFig(figsize = (12,5))

for j in range(0,len(algorithms)):

    pp.plot(gens,fit_mean[:,:,j][0], color = cols[j], label = algorithms[j])
    pp.plot(gens,fit_min[:,:,j][0], color = cols[j], ls = '--', label = None)
    pp.plot(gens,fit_max[:,:,j][0], color = cols[j], ls = ':')
    
legendlist = ['de1220'] + 2* [None] + ['sade'] + 2*[None] + ['gaco'] + \
                2*[None] + ['pso'] + 2*[None] + ['pso_gen'] + 2*[None]
                
legend1 = pp.legend(loc = 1)

pp.legend(['generation mean', 'generation min', 'generation max'], loc=3)
pp.gca().add_artist(legend1)
pp.xlim([1,n_gen-1])
pp.ylim([0,38e3])

ax.set_xlabel('Generation')
ax.set_ylabel('Cost')
pp.tight_layout()
pp.savefig(figure_folder + 'Algocomparison_combined')

#%% avg. similarity plot

fig,ax = PF.MakeFig(figsize=(12,4))

for j in range(0,len(algorithms)):
    pp.plot(gens,resid[j,:], color = cols[j], label = algorithms[j])

pp.xlim([1,n_gen-1])
pp.legend()
ax.set_xlabel('Generation')
ax.set_ylabel('Average cosine similarity')
pp.tight_layout()
pp.savefig(figure_folder + 'Algocomparison_cosine_similarity')
