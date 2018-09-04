"""
Script to compare ASCOT and NUBEAM results
"""
import matplotlib.pyplot as plt
import numpy as np
import transp_output
import ascot_distributions
import scipy.integrate as integrate
import scipy.interpolate as interp

def convert_rhopol_rhotor(locpsi, locq):
    locpsi = locpsi**2
    num=len(locpsi)
    locpsi = interp.interp1d(np.linspace(0,1,len(locpsi)), locpsi)
    locpsi = locpsi(np.linspace(0,1,num))
    locq = interp.interp1d(np.linspace(0,1,len(locq)), locq)
    locq = locq(np.linspace(0,1,num))    
    phi = integrate.cumtrapz(locq,locpsi)
    phi = np.concatenate([[0], phi])
    phi = phi/phi[-1]
    rhophi = phi**0.5
    return rhophi

plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
plt.rc('figure', facecolor='white')
plt.rc('legend', fontsize=20)
nubeam_fname = '/home/vallar/TCV/58823/58823O03.CDF'
ascot_fname = '/home/vallar/ASCOT/runs/TCV/58823/ascot_58823091.h5'
t_ascot = 0.9
d_n = transp_output.output_1d(nubeam_fname)
d_a = ascot_distributions.TCV_1d(ascot_fname)
q = np.flipud(d_a.infile['boozer/qprof'][:])
ind_nubeam = np.argmin(d_n.t[d_n.inj_index]-t_ascot<0)

f_pdep = plt.figure()
ax = f_pdep.add_subplot(111)
ax.plot(d_n.rho, d_n.nb_FKP_vars['pe'][ind_nubeam,:]*1e-3, lw=2.3, color='b', label='Nubeam')

_,[rho,pe], _ = d_a._plot_pe()
rho = convert_rhopol_rhotor(rho, q)
ax.plot(rho, pe, lw=2.3, color='r', label='Ascot')
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'P_e [kW/m^3]')
f_pdep.tight_layout()
plt.show()
