import a5py.ascot5io.ascot5 as a5
from utils.plot_utils import _plot_2d, limit_labels, _plot_1d
import matplotlib.pyplot as plt
import numpy as np

label = ['2D', 'Ripple']
col = ['k', 'r']
f=a5.Ascot('/home/vallar/ASCOT/runs/SA_003/3D/NNB_3D_ripple/ascot.h5')
#3D field
# 0189389985, 100 mrkrs, particle-3505710328, new plasma profiles, wocollisions, gyroorbits, traj
# 1114610632, 30k mrkrs, particle-1080860510, new plasma profiles, wocollisions, gyroorbits
#2D field used is B_2DS-0461682159
# RUN OK f.run_1461926894
# 1315732799, 100 mrkrs,particle-3505710328, new plasma profiles, wocollisions, gyroorbits, traj
# 0374782374, 30k mrkrs, particle-1080860510, new plasma profiles, wocollisions, gyroorbits
wall = f.wall.wall_2D_3851912840.read()
rw = wall['R']; zw=wall['z']

f1 = plt.figure(); ax_rhopitch = f1.add_subplot(111)
f2 = plt.figure(); ax_phitheta = f2.add_subplot(111)
f3 = plt.figure(); ax_etime = f3.add_subplot(111)
f4 = plt.figure(); ax_rz = f4.add_subplot(111)

fig_traj = plt.figure(); axtraj = fig_traj.add_subplot(111); axtraj.set_title(r'Trajectories')
axtraj.plot(rw, zw, 'k', lw=2.3)
#label=['w coll','wo coll']
for ii, run in enumerate([f.run_0374782374, f.run_1114610632]):
    e=run.endstate.read()
    i=run.inistate.read()

    div=np.nanmax(np.abs(e["B_R"]/e["R"] + e["B_R_dR"] + e["B_phi_dphi"]/e["R"] + e["B_z_dz"]))
    print('divergence: ', div)
    ind_wall = np.where(e['endcond']!=10)[0]
    ind_wall = np.where(e['endcond']==8)[0] #wall collision endstate
    #ind_wall = np.where(i['rho']>0.9)[0]
    # B = np.sqrt( np.power(e["B_R"],2) + np.power(e["B_phi"],2) + np.power(e["B_z"],2) )
    # e_perp = e['mu']*B
    # e_para = (e['mass']*1.66054e-27)*0.5*(e['vpar']**2)/1.602e-19
    # e_end = e_perp+e_para
    e_end = run.endstate['energy']/1.602e-19
#    e_end = (e['mass']*1.66054e-27)*0.5*(v_end**2)/1.602e-19
    print('Lost power to wall in case '+str(ii)+' = '+str(np.sum(e_end[ind_wall]*1.602e-19*e['weight'][ind_wall])*1e-6))
    print('Lost markers =', np.size(ind_wall)/np.size(e['endcond']))
    rho = i['rho']  
    v_ini = np.sqrt(i['vR']**2+i['vz']**2+i['vphi']**2)
    e_ini = i['mass']*1.66054e-27*0.5*v_ini*v_ini/1.602e-19
    e_ini = run.inistate['energy']/1.602e-19
    print('P inj in case '+str(ii)+' = '+str(np.sum(i['weight']*e_ini*1.602e-19*1e-6)))
    pitch = i['vpar']/v_ini
    #ax_rhopitch.scatter(rho, pitch, marker='+', color='r')       
    ax_rhopitch.scatter(rho[ind_wall], pitch[ind_wall], color=col[ii])
    # plot trajectories
    try:
        rr = run.orbits.read()
        for index in e['id'][ind_wall[0:5]]:
            _iind = np.where(rr['id']==index)[0]
            axtraj.plot(rr['R'][_iind], rr['z'][_iind])
            #axtraj.plot(rr['R'][_iind][0], rr['z'][_iind][0], 'kx')
        axtraj.plot([i['R'],rr['R'][_iind][0]] , [i['z'],rr['z'][_iind][0]], 'k-')
            
        bb=run.bfield.read()
        _Rb = np.linspace(bb['R_min'], bb['R_max'], bb['n_R'])
        _zb = np.linspace(bb['z_min'], bb['z_max'], bb['n_z'])
        axtraj.contour(_Rb, _zb, bb['psi'], bb['psisepx']) 
        
    except:
        plt.close(fig_traj)
        None


    theta= np.mod(e['pol'], 360.); phi=np.mod(e['phi'], 360)
    _plot_2d(phi[ind_wall],theta[ind_wall], r'$\phi$', r'$\theta$',hist=1)

    time = e['time']
    ax_etime.scatter((e_end[ind_wall]-e_ini[ind_wall])*1e-3, time[ind_wall], label=label[ii])

    R = e['R']; z=e['z']
    #ax_rz.scatter(R[ind_wall], z[ind_wall], color='r')
    ax_rz.plot(rw, zw, 'k', lw=2.3)
    #_plot_2d(rho[ind], pitch[ind], xlabel=r'$\rho$', ylabel=r'$\xi$', scatter=1, ax=ax_rhopitch)
    _plot_2d(i['R'], i['z'], xlabel=r'$R$', ylabel=r'z', scatter=1, ax=ax_rz)    
    #_plot_2d(e_end[ind]*1e-3,e['time'][ind], xlabel=r'$E_{end} [keV]$', ylabel=r't [s]', scatter=1, ax=ax_etime)

    #_plot_1d(e_end[ind]*1e-3, xlabel=r'$E_{end} [keV]$', hist=1, ax=ax_ehisto, color=col[ii])

limit_labels(ax_rhopitch, r'$\rho$', r'$\xi$')
limit_labels(ax_phitheta, r'$\phi$', r'$\theta$')
limit_labels(ax_etime, r'$E_{end}-E_{ini}$[keV]', r'time')
limit_labels(ax_rz, r'$R [m]$', r'z [m]')
ax_etime.legend(loc='best')
#limit_labels(ax_etime, r'E', r't')
f1.tight_layout(); f2.tight_layout(); f3.tight_layout(); f4.tight_layout()
ax_rhopitch.grid('on')
ax_rz.grid('on'); ax_rz.axis('equal')
#ax_etime.grid('on')

plt.show()
