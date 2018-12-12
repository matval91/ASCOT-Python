#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:33:26 2018

@author: vallar
"""
import numpy as np
import utils.plot_utils as ascot_utils

def plot_RZ(data_i, wallrz, surf, xlim, ax=0):
    """plot inistate RZ
    
    Method to plot R vs z of the ionised particles, useful mostly 
    with bbnbi
    
    Parameters:
        | data_i (dict)         : dictionary with inistate. Needs just 'R', 'z' keys
        | wallrz (array)        : 2D array with [R,Z] of the wall. Will plot RZ wall section
        | surf (array)          : 2D array with [R,Z, surf] of the flux surfaces. 
                                    Will plot flux surfaces wrt R,Z 
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)

    Arguments:
        None    
    """
    # bypass problem with Rprt (sometimes it is only 999)
    try:
        x=data_i['Rprt']
        y=data_i['zprt']
    except:
        x=data_i['R']
        y=data_i['z']

    if np.mean(x)==999. or np.mean(x)==-999.:
        x=data_i['R']
        y=data_i['z']
        
    xlab = r'R [m]'; ylab = r'z [m]'
    
    ascot_utils._plot_2d(x, y, xlabel=xlab, ylabel=ylab, title='RZ ionization',\
             wallrz=wallrz, surf=surf, ax=ax, xlim=xlim, scatter=1)


#    if 'bbnbi' in self._fname and shpart!=0:
#        ax = plt.gca()
#        R = data_e['R']
#        ind = np.where(R<4.5)
#        R=R[ind]
#        phi = data_e['phi'][ind]
#        x,y = R*np.cos(phi), R*np.sin(phi)
#        z = data_e['z'][ind]
#        _plot_2d(R,z, scatter=1, ax=ax)
    
    
def plot_XY(data_i, wallxy, R0=0, ax=0, shpart=0):
    """plot inistate RZ
    
    Method to plot XY of ionisation, without difference between the beams
    
    Parameters:
        | data_i (dict)         : dictionary with inistate. Needs just 'R', 'z' keys
        | wallxy (array)        : 2D array with [R,Z] of the wall. Will plot circles with
                                    min(R) and max(R) (default 0)
        | R0 (int)              : axis R0
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)

    Arguments:
        None    
    """

    try:
        R=data_i['Rprt']
        phi = data_i['phiprt']
    except:
        R=data_i['R']
        phi = data_i['phi']

    if np.abs(np.mean(R))==999.:
        R=data_i['R']
        phi = data_i['phi']

    x = R*np.cos(phi)
    y = R*np.sin(phi)
    xlab = r'X [m]';    ylab = r'Y [m]' 
    
    ascot_utils._plot_2d(x, y, xlabel=xlab, ylabel=ylab, title='XY Ionization',\
             wallxy=wallxy, R0=R0, ax=ax)