import a5py.ascot5io.ascot5 as a5
from utils.plot_utils import _plot_2d, limit_labels, _plot_1d
import matplotlib.pyplot as plt
import numpy as np

label = ['2D', 'Ripple']
col = ['k', 'r']
folder='/home/vallar/ASCOT/runs/SA_003/3D/NNB_3D_ripple/'
#folder += 'high_resolution/diffusive/'
f=a5.Ascot(folder+'ascot.h5')

#3D field used is B_3DS-3439715758
# RUN OK 0047991044
#run with high resolution for orbits: f.run_1553361412

f1 = plt.figure(); ax_rhopitch = f1.add_subplot(111)
f2 = plt.figure(); ax_Rz = f2.add_subplot(111)
f3 = plt.figure(); ax_ehisto = f3.add_subplot(111)
w=f.wall.wall_2D_3851912840.read()
Rw=w['R']; zw=w['z']


for ii, run in enumerate([f.run_1620750822]):#0291578081]):#, f.run_0047991044]):
    e=run.endstate.read()
    endcond=e['endcond']
    ind=np.where(e['endcond']!=10)[0] #wall collision endstate
    R=e['R']; z=e['z']
    ax_Rz.scatter(R[ind],z[ind], color=col[ii])
    ax_Rz.plot(Rw, zw, 'k')
    e['vz'] = e['vz']*1.66054e-27
    v_end = np.sqrt(e['vR']**2+e['vz']**2+e['vphi']**2)
    e_end = e['mass']*1.66054e-27*0.5*v_end*v_end/1.602e-19
    print('Lost power in case '+str(ii)+' = '+str(np.sum(e_end[ind]*1.602e-19*e['weight'][ind])*1e-6))
    i=run.inistate.read()
    print('P inj in case '+str(ii)+' = '+str(np.sum(i['weight']*500e3*1.602e-19*1e-6)))

    rho = i['rho']
    R = i['R']; z=i['z']
    v_ini = np.sqrt(i['vR']**2+i['vz']**2+i['vphi']**2)
    pitch = i['vpar']/v_ini
    B = np.sqrt(i['B_phi']**2+i['B_R']**2+i['B_z']**2)
    ax_rhopitch.scatter(rho[ind], pitch[ind], color=col[ii])
    #_plot_2d(rho[ind], pitch[ind], xlabel=r'$\rho$', ylabel=r'$\xi$', scatter=1, ax=ax_rhopitch)
    _plot_2d(R[ind], z[ind], xlabel=r'$R$', ylabel=r'B', scatter=1, ax=ax_Rz)    
    #_plot_2d(e_end[ind]*1e-3,e['time'][ind], xlabel=r'$E_{end} [keV]$', ylabel=r't [s]', scatter=1, ax=ax_etime)

    #_plot_1d(e_end[ind]*1e-3, xlabel=r'$E_{end} [keV]$', hist=1, ax=ax_ehisto, color=col[ii])

limit_labels(ax_rhopitch, r'$\rho$', r'$\xi$')
#limit_labels(ax_etime, r'E', r't')
f1.tight_layout()
f2.tight_layout()
ax_rhopitch.grid('on')
#ax_etime.grid('on')


plt.show()
