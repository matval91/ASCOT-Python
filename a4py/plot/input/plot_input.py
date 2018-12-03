#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
These functions should be used only to plot input to ascot4


Created on Wed Nov 28 18:51:00 2018

@author: vallar
"""
import utils.ascot_utils as ascot_utils
import matplotlib.pyplot as plt

colours = ['k','b','r','c','g']

def plot_profiles(prof, fig=0, title=''):
    """Plot the profiles
    
    This function makes a plot with ne, Te, Ti, ni(eventually nimp) on 4 different frames

    Parameters:
        | prof (object): object created using ascot_prof
        |  f (object): the plt.figure object where to plot (useful to overplot). Undefined is initialized to 0 and it means to do a new figure
        |  title (str): title of the figure
    Return:
        None

    """
    ascot_utils.common_style()
    

    #=====================================================================================
    if fig==0:
        w, h = plt.figaspect(0.8)
        fig=plt.figure(title, figsize=(10,8))
        axte = fig.add_subplot(221)
        axne = fig.add_subplot(222, sharex=axte)
        axti = fig.add_subplot(223, sharey=axte)
        axni = fig.add_subplot(224, sharey=axne)
    else:
        axte = fig.axes[0]
        axne = fig.axes[1]
        axti = fig.axes[2]
        axni = fig.axes[3]
    #axvt = fig.add_subplot(325)
    lw=4
    axte.plot(prof.rho, prof.te*1e-3,'k', linewidth=lw)
    axne.plot(prof.rho, prof.ne,'k', linewidth=lw)
    axti.plot(prof.rho, prof.ti*1e-3,'k', linewidth=lw)
    #axvt.plot(self.rho, self.vt,'k', linewidth=lw)
    for i in range(prof.nion):
        if prof.A[i]==2.:
            label='D'
        elif prof.A[i]==12:
            label='C'
        else:
            label='UNSPEC'
        axni.plot(prof.rho, prof.ni[i,:],colours[i], \
                  label=label, linewidth=lw)
    axni.legend(loc='best')


    if fig == 0:
        ascot_utils.limit_labels(axte, r'$\rho$', r'$T_e$ [keV]')
        ascot_utils.limit_labels(axte, r'$\rho$', r'$n_e$ [1/$m^3$]')
        ascot_utils.limit_labels(axte, r'$\rho$', r'$T_i$ [keV]')
        ascot_utils.limit_labels(axte, r'$\rho$', r'$n_i$ [1/$m^3$]')

        fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    plt.show()