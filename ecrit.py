#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:51:24 2018

@author: vallar
"""
import h5py
import numpy as np
import scipy.interpolate as interpolate

def ecrit(inf_name='ascot.h5', E0=500000):
    """ecrit calculation
    Calculates critical energy profiles
    Ec = 14.8*te*(A**(1.5)/ne*summ)**(2./3.)
    ts = 6.28e14*(A*te^1.5)/(Z^2*ne*lnlambda)
    note: tega
    """
    print("Calculating ts with E0="+str(E0)) 
    infile=h5py.File(inf_name)
    rho = infile['plasma/1d/rho'][:]
    volumes = infile['distributions/rhoDist/shellVolume'].value
    _volumes = interpolate.interp1d(np.linspace(0,1,len(volumes)), volumes)
    volumes=_volumes(rho[rho<1])
    ind=rho<1;
    te = infile['plasma/1d/te'][ind]
    ne = infile['plasma/1d/ne'][ind]
    Ai = infile['plasma/anum'][:]
    Zi = infile['plasma/znum'][:]
    nimp = infile['plasma/1d/ni'][ind,:]
    rho = rho[ind]
    A = infile['species/testParticle/anum'][0]
    Z = infile['species/testParticle/znum'][0]
    summ = np.sum(np.multiply(nimp, Zi**2/Ai), axis=1)
    Ec = 14.8*te*(A**(1.5)/ne*summ)**(2./3.)
    param_ec = interpolate.interp1d(rho,Ec)
    ec_mean = np.trapz(param_ec(rho)*volumes)/np.sum(volumes)
    
    #Spitzer slowing-down time
    ts = 6.28e14*A*te**1.5/(Z**2*ne*17.)
    taus=ts/3.*np.log((1+(E0/Ec)**1.5))
    param_taus = interpolate.interp1d(rho,taus)
    taus_mean = np.trapz(param_taus(rho)*volumes)/np.sum(volumes)
    return param_ec, ec_mean, param_taus, taus_mean
