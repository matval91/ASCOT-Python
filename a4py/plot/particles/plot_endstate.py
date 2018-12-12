#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:25:45 2018

@author: vallar
"""
import numpy as np
import matplotlib.pyplot as plt
import utils.plot_utils as ascot_utils

def plot_histo_wall(data_e, R0=0, z0=0, ax=0, lastpoint=1):
    """ plot wall losses
    
    Histogram of deposition to the wall
    
    Parameters:
        | data_e (dict)         : dictionary with endstate. Needs just 'endcond', 
                                    'R', 'z' keys
        | R0 (int)              : axis R0, needed for theta
        | z0 (int)              : axis z0, needed for theta
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)
        | lastpoint (str)       : if 0, doesn't plot last point in colorbar (default=1)

    Arguments:
        None
    """
    ascot_utils.common_style()
    
    # filtering on data
    ind = np.where(data_e['endcond']== 3)[0] #wall
    r = data_e['R'][ind]
    z = data_e['z'][ind]

    theta = np.arctan2(z-z0,r-R0)
    phi = data_e['phi'][ind]

    # defining labels
    angles_ticks=[-3.14, -1.57, 0., 1.57, 3.14]
    angles_labels=[r'-$\pi$',r'-$\pi/2$',r'0',r'$\pi/2$',r'$\pi$']
    
    ascot_utils._plot_2d(phi, theta, xlabel=r'$\phi$ [rad]',ylabel=r'$\theta$ [rad]', \
             hist=1, xlim=[-3.14, 3.14],ylim=[-3.14, 3.14], cblabel='# markers', \
             lastpoint=lastpoint)
    # Setting labels to radians values
    ax=plt.gca()
    ax.set_xticks(angles_ticks); ax.set_xticklabels(angles_labels)
    ax.set_yticks(angles_ticks); ax.set_yticklabels(angles_labels)

    # Now plotting theta distribution
    ascot_utils._plot_1d(theta, xlabel=r'$\theta$ [rad]', hist=1)
    ax=plt.gca()
    ax.set_xticks(angles_ticks); ax.set_xticklabels(angles_labels)
    
    # Now plotting phi distribution
    ascot_utils._plot_1d(phi, xlabel=r'$\phi$ [rad]', hist=1)
    ax=plt.gca()
    ax.set_xticks(angles_ticks); ax.set_xticklabels(angles_labels)
    
def plot_rhopitch(data_e, data_i):
    """ plot rhopitch of losses

    Method to plot rho, pitch histo for losses
    
    Parameters:
        | data_e (dict)         : dictionary with endstate. Needs 'endcond' keys
        | data_e (dict)         : dictionary with inistate. Needs 'R', 'z' keys
        
    Returns:
        None
    """
    ind = np.where(data_e['endcond']==3)[0]
    rho = data_i['rho'][ind]
    pitch = data_i['pitch'][ind]
    xlabel = r'$\rho$'; ylabel=r'$\xi=v_\parallel/v$'
    ascot_utils._plot_2d(rho, pitch, hist=1, ylim=[-1, 1], xlabel=xlabel, ylabel=ylabel)