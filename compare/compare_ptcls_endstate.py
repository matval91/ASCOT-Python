# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 09:53:29 2018

Comparing two simulations, in this case for particles end status
Created to compare the two results with r9401 and r9512

@author: vallar
"""

import ascot_particles
import matplotlib.pyplot as plt
import numpy as np
import math

col=['k','r','b']

SA005dir = '/home/vallar/ASCOT/runs/JT60SA/005/'
p08=ascot_particles.h5_particles(SA005dir+'ascot_005008.h5') #run with r9215
p22=ascot_particles.h5_particles(SA005dir+'ascot_005022.h5') #run with r9401
p23=ascot_particles.h5_particles(SA005dir+'ascot_005023.h5') #run with r9513

f = plt.figure()
axpitch  = f.add_subplot(221)
axpitch.set_title('Pitch')
axenergy = f.add_subplot(222)
axenergy.set_title('Energy')
axxy = f.add_subplot(223)
axxy.set_title('XY')

for ii, obj in enumerate([p08, p22, p23]):
    ind = np.where(obj.data_e['endcond']== 4)[0] #thermalized

    pitchi = obj.data_i['pitch'][ind]
    energyi = obj.data_i['energy'][ind]
    pitch = obj.data_e['pitch'][ind]
    energy = obj.data_e['energy'][ind]
    vr = obj.data_e['vR'][ind]
    vz = obj.data_e['vz'][ind]
    vphi = obj.data_e['vphi'][ind]
    norm_v = np.sqrt(vr**2+vz**2+vphi**2)
    r = obj.data_e['R'][ind]
    z = obj.data_e['z'][ind]
    theta = np.arctan2(z,r-3.)*180./math.pi
    phi = obj.data_e['phi'][ind]
    x = r*np.cos(phi); y=r*np.sin(phi) 
    
    axpitch. hist(pitch , normed='false', alpha=0.5, color=col[ii])
    axenergy.hist(energy, normed='false', alpha=0.5, color=col[ii])

plt.show()