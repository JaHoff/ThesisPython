# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 16:52:07 2020

@author: USER
"""

from os import getcwd
from os import chdir
import sys 
import numpy as np
import matplotlib.pyplot as pp
# from adjustText import adjust_text

sys.path.append('../')
chdir('../')
import DataImport as DI
import GeneralFunctions as GF
import LazyLib as LL
import PlottingFunctions as PF
chdir('PerturbationAnalysis')

figfolder = LL.folder_dict["NumericalSim"]

def makeHodograph(data, planet = "", ax1 = 'x', ax2 = 'y'):
    
    indexes = {'x':0 ,
               'y':1 ,
               'z':2 }
    
    
    angles = np.arctan(data[:,indexes[ax2]]/data[:,indexes[ax1]]) 
    # * np.sign(data[:,indexes[ax2]]) * np.sign(data[:,indexes[ax1]])
    # angles = np.abs(angles)
    angles[(data[:,indexes[ax1]] < 0) * (data[:,indexes[ax2]]>0)] += np.pi
    angles[(data[:,indexes[ax1]] < 0) * (data[:,indexes[ax2]]<0)] += np.pi
    angles[(data[:,indexes[ax1]] > 0) * (data[:,indexes[ax2]]<0)] += 2*np.pi
    # if ax1 == 'y' and ax2 == 'z':
    #     angles = np.arctan(data[:,indexes[ax2]]/data[:,indexes[ax1]])
    # else:
    #     angles = np.arccos(data[:,indexes[ax1]])
        
    # angles[data[:,indexes[ax2]]<0] = 2*np.pi-angles[data[:,indexes[ax2]]<0] 
    
    boxes = np.arange(0,2*np.pi,2*np.pi/10)
    
    vals = np.histogram(angles,120,(0,2*np.pi))
    
    maxval = np.max(vals[0])
    meanval = np.mean(vals[0][vals[0] > 0])
    percentval = 100/np.sum(vals[0])
    zeros = np.empty_like(data[:,indexes[ax1]])*0
    
    pp.figure()
    
    for i in range(0,len(vals[0])):
        ind = np.mod(i+1,len(vals[0]))
        x = np.array([0, vals[0][i]*np.cos(vals[1][i]), vals[0][ind]*np.cos(vals[1][i+1]), 0])*percentval
        y = np.array([0, vals[0][i]*np.sin(vals[1][i]), vals[0][ind]*np.sin(vals[1][i+1]), 0])*percentval
        pp.plot(x,y, c= 'k')
    
    
    pp.title(f"Relative pull strength: {planet}")
    pp.xlabel(f"Relative pull strength, barycentric {ax1}")
    pp.ylabel(f"Relative strength, barycentric {ax2}")
    pp.xlim([-1.1, 1.1])
    pp.ylim([-1.1, 1.1])
    
    pp.savefig(figfolder + f"relative_pull_hodograph_{planet}_{ax1}_{ax2}.png")
    return vals,angles
    
def gravitation(rel_pos, m1, m2= 5):
    G = LL.Constants["G"]
    R = np.linalg.norm(rel_pos,axis=1)
    F = G*m1*m2/(R**2)
    return F
    
#%% 
chdir("../")
cwd = getcwd()
append = cwd + "/Data/"

print("Start synchronizing data files to folder %s" % append)
LL.SyncDataFiles(targetDir = append)
print("Synchronized data files")

# x,y,z,vx,vy,vz, R, V = (np.ndarray([]) for i in range(0,8))

datafolder = cwd + '/Data/DefaultModel/'

#%%
PF.CloseAll()
pos = DI.ImportDataFile(datafolder + 'relative_body_pos_365d.dat', delim = ',')

moon = pos[:,:6]
sun = pos[:,6:9]
mercury = pos[:,9:12]
mars = pos[:,12:15]
venus = pos[:,15:18]
jupiter = pos[:,18:21]
saturn =  pos[:,21:]


b_sun = DI.ConvertToBarycentricFrame(moon,sun, exportL4 = True, fixedL4 = True)
b_mercury = DI.ConvertToBarycentricFrame(moon,mercury)
b_mars = DI.ConvertToBarycentricFrame(moon,mars)
b_venus = DI.ConvertToBarycentricFrame(moon,venus)
b_jupiter = DI.ConvertToBarycentricFrame(moon,jupiter)
b_saturn = DI.ConvertToBarycentricFrame(moon,saturn)

l4 = b_sun[3]

b_sun = np.concatenate((b_sun[0], b_sun[1],b_sun[2]), axis=1)
b_mercury = np.concatenate((b_mercury[0], b_mercury[1],b_mercury[2]), axis=1)
b_mars = np.concatenate((b_mars[0], b_mars[1],b_mars[2]), axis=1)
b_venus = np.concatenate((b_venus[0], b_venus[1],b_venus[2]), axis=1)
b_jupiter = np.concatenate((b_jupiter[0], b_jupiter[1],b_jupiter[2]), axis=1)
b_saturn = np.concatenate((b_saturn[0], b_saturn[1],b_saturn[2]), axis=1)

#%% relative positions to L4

rel_sun = b_sun - l4
rel_mer = b_mercury - l4
rel_mar = b_mars - l4
rel_ven = b_venus - l4
rel_jup = b_jupiter - l4
rel_sat = b_saturn - l4

rel_sun /= np.linalg.norm(rel_sun,axis = 1).reshape(17533,1)
rel_mer /= np.linalg.norm(rel_mer,axis = 1).reshape(17533,1)
rel_mar /= np.linalg.norm(rel_mar,axis = 1).reshape(17533,1)
rel_ven /= np.linalg.norm(rel_ven,axis = 1).reshape(17533,1)
rel_jup /= np.linalg.norm(rel_jup,axis = 1).reshape(17533,1)
rel_sat /= np.linalg.norm(rel_sat,axis = 1).reshape(17533,1)

rel_ear = -l4 / np.linalg.norm(l4,axis=1).reshape(17533,1)

#%%

makeHodograph(rel_sun, "Sun", 'x', 'y')
makeHodograph(rel_mer, "Mercury",'x', 'y')
makeHodograph(rel_mar, "Mars",'x', 'y')
makeHodograph(rel_ven, "Venus",'x', 'y')
makeHodograph(rel_jup, "Jupiter",'x', 'y')
makeHodograph(rel_sat, "Saturn", 'x', 'y')


makeHodograph(rel_sun, "Sun", 'y', 'z')
makeHodograph(rel_mer, "Mercury",'y', 'z')
makeHodograph(rel_mar, "Mars",'y', 'z')
makeHodograph(rel_ven, "Venus",'y', 'z')
makeHodograph(rel_jup, "Jupiter",'y', 'z')
makeHodograph(rel_sat, "Saturn", 'y', 'z')
# pp.quiver(zeros,zeros,np.cos(vals[1][:-1]),np.sin(vals[1][:-1]), scale = vals[0])
# pp.xlim([-400,400])
# pp.ylim([-400,400])

#%% Make a plot of planet positions

AU = LL.Constants["AU"]
fig = pp.figure(figsize=(12,3.3))
ax = fig.gca()


pp.scatter(0 , 0, color='y', marker = '*')
pp.scatter(-sun[0,1]/AU, -sun[0,0]/AU,  color = 'b')
# pp.scatter( (moon[0,1]-sun[0,1])/AU, (moon[0,0]-sun[0,0])/AU, color = 'b', s = 3)
pp.scatter((mercury[0,1]-sun[0,1])/AU,(mercury[0,0]-sun[0,0])/AU, color = 'k')
pp.scatter( (venus[0,1]-sun[0,1])/AU,(venus[0,0]-sun[0,0])/AU, color = 'g')
# pp.scatter((mars[0,0]-sun[0,0])/AU, (mars[0,1]-sun[0,1])/AU, color='r')
# pp.scatter((jupiter[0,0]-sun[0,0])/AU, (jupiter[0,1]-sun[0,1])/AU, color = 'y')
# pp.scatter((saturn[0,0]-sun[0,0])/AU, (saturn[0,1]-sun[0,1])/AU, color = 'y')

r = np.sqrt(sun[0,0]**2 + sun[0,1]**2)/AU
ax.add_artist(pp.Circle((0,0),r,color = 'b',fill = False))
r = np.sqrt((mercury[0,0]-sun[0,0])**2 + (mercury[0,1]-sun[0,1])**2)/AU
ax.add_artist(pp.Circle((0,0),r,color = 'k',fill = False))
r = np.sqrt((venus[0,0]-sun[0,0])**2 + (venus[0,1]-sun[0,1])**2)/AU
ax.add_artist(pp.Circle((0,0),r,color = 'g',fill = False))
# r = np.sqrt((mars[0,0]-sun[0,0])**2 + (mars[0,1]-sun[0,1])**2)/AU
# ax.add_artist(pp.Circle((0,0),r,color = 'r',fill = False))
# r = np.sqrt((jupiter[0,0]-sun[0,0])**2 + (jupiter[0,1]-sun[0,1])**2)/AU
# ax.add_artist(pp.Circle((0,0),r,color = 'y',fill = False))
# r = np.sqrt((saturn[0,0]-sun[0,0])**2 + (saturn[0,1]-sun[0,1])**2)/AU
# ax.add_artist(pp.Circle((0,0),r,color = 'y',fill = False))

pp.xlim([-1, 0.1])
pp.ylim([-0.25, .25])
pp.xlabel('J2000-y [AU]')
pp.ylabel('J2000-x [AU]')
pp.legend(['Sun', 'Earth', 'Mercury', 'Venus' ])
pp.tight_layout()
pp.savefig(figfolder + 'Initial_solar_system')
#%%  Co-linearity plot
x = np.linspace(0,365,len(sun[:,0]))
fig,ax = PF.MakeFig(figsize=(12,4))
pp.plot(x, np.linalg.norm(b_mercury-l4,axis=1)/AU)
pp.plot(x, np.linalg.norm(b_venus-l4,axis=1)/AU)
pp.plot(x, np.linalg.norm(b_mars-l4,axis=1)/AU)
pp.plot(x, np.linalg.norm(b_jupiter-l4,axis=1)/AU)
pp.plot(x, np.linalg.norm(b_saturn-l4,axis=1)/AU)

pp.xlabel('Time in orbit [days]')
pp.ylabel('Distance from L4 [AU]')

#%%  Gravitational pull plot

m_merc =3.3011e23
m_venus = 4.8675e24
m_mars = 6.4171e23
m_jup = 1.8982e27
m_sat = 5.6834e26
m_sun = 1.9884e30

fig,ax = PF.MakeFig(figsize=(12,3.3))
pp.plot(x, gravitation(b_mercury-l4, m_merc)*1e3)
pp.plot(x, gravitation(b_venus-l4, m_venus)*1e3)
pp.plot(x, gravitation(b_mars-l4, m_mars)*1e3)
pp.plot(x, gravitation(b_jupiter-l4, m_jup)*1e3)
pp.plot(x, gravitation(b_saturn-l4, m_sat)*1e3, ls = '--')
pp.plot(x, gravitation(b_sun-l4, m_sun)*1e3, ls='--')


pp.xlabel('Time in orbit [days]')
pp.ylabel('Point mass gravitation [mN]')
pp.legend(['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Sun'])
pp.xlim([0,365])
pp.yscale('log')

pp.tight_layout()
pp.savefig(figfolder + 'grav_pull')

#%% Orbit resonance investigation



datafolder = cwd + '/Data/PerturbationAnalysis/'
filename = datafolder + "propagationHistory_" + str(0) + '_100km_nomoon' + ".dat"
t, xd, yd, zd, vxd, vyd, vzd,mn = DI.ImportPropagationHistory(filename,0, False)
res_bor = 0.91

y_orb = (xd[:,0]-l4[:,0])/np.max(np.abs(xd[:,0]-l4[:,0]))
y_ven = rel_ven[:,0]/np.max(np.abs(rel_ven[:,0]))
y_jup = rel_jup[:,0]/np.max(np.abs(rel_jup[:,0]))

res_ven = np.empty_like(y_orb)*np.nan
# res_ven[y_orb*y_ven > res_bor] = y_orb[y_orb*y_ven > res_bor]
# res_jup = np.empty_like(y_orb)*np.nan
# res_jup[y_orb*y_jup > res_bor] = y_orb[y_orb*y_jup > res_bor]

fig,ax = pp.subplots(3,1,figsize=(12,9),sharex=True)

ax[0].plot(x,y_orb)
ax[0].plot(x,y_ven, color = 'r', ls = '--')
ax[0].plot(x,y_jup, color = 'g', ls = '--')
# ax[0].scatter(x,res_ven,marker='d', zorder = 20, color = 'r')
# ax[0].scatter(x,res_jup,marker='d', zorder = 20, color = 'g')
ax[0].set_ylabel('Normalized relative x')
# pp.title('x-resonance')
# pp.legend(['x-orbit', 'Venus - rel. x', 'Jupiter - rel x'])

y_orb = (yd[:,0]-l4[:,1])/np.max(np.abs(yd[:,0]-l4[:,1]))
y_ven = rel_ven[:,1]/np.max(np.abs(rel_ven[:,1]))
y_jup = rel_jup[:,1]/np.max(np.abs(rel_jup[:,1]))

res_ven = np.empty_like(y_orb)*np.nan
# res_ven[y_orb*y_ven > res_bor] = y_orb[y_orb*y_ven > res_bor]
# res_jup = np.empty_like(y_orb)*np.nan
# res_jup[y_orb*y_jup > res_bor] = y_orb[y_orb*y_jup > res_bor]

ax[1].plot(x,y_orb)
ax[1].plot(x,y_ven, color = 'r',ls = '--')
ax[1].plot(x,y_jup, color = 'g', ls = '--')
# ax[1].scatter(x,res_ven,marker='d', zorder = 20, color = 'r')
# ax[1].scatter(x,res_jup,marker='d', zorder = 20, color = 'g')
ax[1].set_ylabel('Normalized relative y')
ax[1].set_xlim([0,365])


y_orb = (zd[:,0]-l4[:,2])/np.max(np.abs(zd[:,0]-l4[:,2]))
y_ven = rel_ven[:,2]/np.max(np.abs(rel_ven[:,2]))
y_jup = rel_jup[:,2]/np.max(np.abs(rel_jup[:,2]))

res_ven = np.empty_like(y_orb)*np.nan
res_ven[y_orb*y_ven > res_bor] = y_orb[y_orb*y_ven > res_bor]
res_jup = np.empty_like(y_orb)*np.nan
res_jup[y_orb*y_jup > res_bor] = y_orb[y_orb*y_jup > res_bor]


ax[2].plot(x,y_orb)
ax[2].plot(x,y_ven, color = 'r',ls = '--')
ax[2].plot(x,y_jup, color = 'g',ls = '--')
# ax[2].scatter(x,res_ven,marker='d', zorder = 20, color = 'r')
# ax[2].scatter(x,res_jup,marker='d', zorder = 20, color = 'g')
ax[2].set_xlim([0,365])
ax[2].set_ylabel('Normalized relative z')
ax[2].set_xlabel('Time [days]')

fig.legend(['Relative orbit $S_1$', 'Venus - relative', 'Jupiter - relative'],loc=(0.83,0.07))
pp.tight_layout()

pp.savefig(figfolder + "resonance_plot")