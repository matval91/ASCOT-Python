#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
These functions should be used only to plot input to ascot4


Created on Wed Nov 28 18:51:00 2018

@author: vallar
"""
import utils.plot_utils as au
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
    au.common_style()
    overplot=False
    ls = '-'
    if fig==0:
        w, h = plt.figaspect(0.8)
        fig=plt.figure(title, figsize=(10,8))
        axte = fig.add_subplot(221)
        axne = fig.add_subplot(222, sharex=axte)
        axti = fig.add_subplot(223, sharey=axte, sharex=axte)
        axni = fig.add_subplot(224, sharey=axne, sharex=axte)
    else:
        overplot=True
        ls='--'
        axte = fig.axes[0]
        axne = fig.axes[1]
        axti = fig.axes[2]
        axni = fig.axes[3]
    #axvt = fig.add_subplot(325)
    lw=2.3
    axte.plot(prof.rho, prof.te*1e-3,'k', linewidth=lw, linestyle=ls)
    axne.plot(prof.rho, prof.ne,'k', linewidth=lw, linestyle=ls)
    axti.plot(prof.rho, prof.ti*1e-3,'k', linewidth=lw, linestyle=ls)
    #axvt.plot(self.rho, self.vt,'k', linewidth=lw)
    for i in range(prof.nion):
        if prof.A[i]==2.:
            label='D'
        elif prof.A[i]==12:
            label='C'
        else:
            label='UNSPEC'
        axni.plot(prof.rho, prof.ni[i,:],colours[i], \
                  label=label, linewidth=lw, linestyle=ls)
    axni.legend(loc='best')
    axte.set_xlim([0,1.])
    if ~overplot:
        au.limit_labels(axte, r'$\rho_{pol}$', r'$T_e$ [keV]')
        au.limit_labels(axne, r'$\rho_{pol}$', r'$n_e$ [1/$m^3$]')
        au.limit_labels(axti, r'$\rho_{pol}$', r'$T_i$ [keV]')
        au.limit_labels(axni, r'$\rho_{pol}$', r'$n_i$ [1/$m^3$]')

        fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    plt.show()

def plot_Bfield(B, fig=0, title=''):
    """plot of 2D psi, q, bphi

    Method to plot the values (2D psi, q, bphi) and check the magnetic field we are looking at
    
    Parameters:
        f (obj): figure object where to plot. if undefined, fig=0
    Attributes:
        None

    """
    au.common_style()
    
    if fig==0:
        fig = plt.figure(figsize=(20, 8))
        fig.text(0.01, 0.01, '')
        
    ax2d = fig.add_subplot(131)

    CS = ax2d.contour(B.R,B.z, B.psi, 30)
    plt.contour(B.R,B.z, B.psi, [B.psiedge], colors='k', linewidths=3.)
    plt.colorbar(CS)
    au.limit_labels(ax2d, r'R [m]', r'z [m]', M=3)

    if B.R_w[0]!=0:
        ax2d.plot(B.R_w, B.z_w, 'k',linewidth=2)
    ax2d.axis('equal')
    
    axq = fig.add_subplot(132)
    axq.plot(B.rhopsi, B.q, lw=2.3, color='k')
    au.limit_labels(axq, r'$\rho_{POL}$', r'q', M=3)

    axf = fig.add_subplot(133)
    #axf.plot(B.R_eqd, B.eqdsk.T)
    if len(B.Bphi.shape)>1:
        axf.plot(B.R, B.Bphi[int(len(B.z)/2.),:], lw=2.3, color='k')
    else:
        axf.plot(B.R, B.Bphi, lw=2.3, color='k')
    au.limit_labels(axf, r'R [m]', r'$B_\phi$ [T]', M=3)

    fig.tight_layout()
    
    plt.show()
