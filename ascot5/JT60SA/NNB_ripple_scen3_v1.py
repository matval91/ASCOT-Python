import a5py.ascot5io.ascot5 as a5
from utils.plot_utils import _plot_2d, limit_labels, _plot_1d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


label = ['2D', 'Ripple']
col = ['k', 'r']
f=a5.Ascot('/home/vallar/ASCOT/runs/SA_003/3D/NNB_3D_ripple/ascot.h5')
wall = f.wall.wall_2D_3851912840.read()
rw = wall['R']; zw=wall['z']
#3D field
# 0189389985, 100 mrkrs, particle-3505710328, new plasma profiles, wocollisions, gyroorbits, traj
# 1114610632, 30k mrkrs, particle-1080860510, new plasma profiles, wocollisions, gyroorbits

#2D field used is B_2DS-0461682159
# RUN OK f.run_1461926894
# 1315732799, 100 mrkrs,particle-3505710328, new plasma profiles, wocollisions, gyroorbits, traj
# 0374782374, 30k mrkrs, particle-1080860510, new plasma profiles, wocollisions, gyroorbits

runs=[f.run_1114610632, f.run_0374782374]

def plot_rz(i, e, ind_wall):
    """
    Plot of initial state r,z of particles, considering also those hitting the wall
    """
    f=plt.figure(figsize=(5,8)); ax_rz=f.add_subplot(111)
    _plot_2d(i['R'], i['z'], xlabel=r'$R$', ylabel=r'z', scatter='b', ax=ax_rz)
    _plot_2d(i['R'][ind_wall], i['z'][ind_wall], scatter='r', ax=ax_rz)
    #_plot_2d(e['R'][ind_wall], e['z'][ind_wall], scatter='g', ax=ax_rz)
    #patch = [mpatches.Patch(color='blue', label='Confined'),\
    #         mpatches.Patch(color='red', label='Lost'),\
    #         mpatches.Patch(color='red', label='Wall')]    
    patch = [mpatches.Patch(color='blue', label='Confined'),\
             mpatches.Patch(color='red', label='Lost')]
    ax_rz.plot(rw, zw, 'k', lw=2.3)
    limit_labels(ax_rz, r'$R [m]$', r'z [m]')
    ax_rz.axis('equal'); f.tight_layout()
    plt.legend(handles=patch)

def plot_rhopitch(i, ind_wall):
    """
    Plot of initial state rho,pitch of particles, considering also those hitting the wall
    """
    f=plt.figure(); ax_rhopitch=f.add_subplot(111)
    _plot_2d(i['rho'], i['pitch'], xlabel=r'$\rho$', ylabel=r'$\xi$', scatter='b', ax=ax_rhopitch)
    _plot_2d(i['rho'][ind_wall], i['pitch'][ind_wall], scatter='r', ax=ax_rhopitch)
    limit_labels(ax_rhopitch, r'$\rho$', r'$\xi$')
    patch = [mpatches.Patch(color='blue', label='Confined'),\
    mpatches.Patch(color='red', label='Lost')]

    plt.legend(handles=patch, loc='best')
    f.tight_layout()

def plot_phitheta(e, ind_wall):
    """
    """
    f = plt.figure(); ax_phitheta = f.add_subplot(111)
    theta= np.arctan2(e['z'],e['R']-2.97)*180./np.pi; phi=np.mod(e['phi'], 360)
    _plot_2d(phi[ind_wall],theta[ind_wall], r'$\phi$', r'$\theta$',hist=1, ax=ax_phitheta)
    limit_labels(ax_phitheta, r'$\phi$', r'$\theta$')
    f.tight_layout()

def plot_phi(e, ind_wall):
    """
    """
    f = plt.figure(); ax_phitheta = f.add_subplot(111)
    phi=np.mod(e['phi'], 360)
    _plot_1d(phi[ind_wall], r'$\phi$',hist=1, ax=ax_phitheta)
    limit_labels(ax_phitheta, r'$\phi$', r'Markers (AU)')
    f.tight_layout()

def plot_E(e, ind_wall):
    """
    """
    f = plt.figure(); ax_phitheta = f.add_subplot(111)
    _plot_1d(e[ind_wall]*1e-3/1.602e-19, r'E [keV]',hist=1, ax=ax_phitheta)
    limit_labels(ax_phitheta, r'E [keV]', r'Markers (AU)')
    f.tight_layout()

def plot_trajectory(rr,e,ind_wall, run):
    """
    """
    fig_traj = plt.figure(figsize=(5,8)); axtraj = fig_traj.add_subplot(111); axtraj.set_title(r'Trajectories')
    fig_traj2 = plt.figure(); axtraj2 = fig_traj2.add_subplot(111); axtraj2.set_title(r'Trajectories')
    axtraj.plot(rw, zw, 'k', lw=2.3)
    for index in e['id'][ind_wall[0:3]]:
        _iind = np.where(rr['id']==index)[0]
        axtraj.plot(rr['R'][_iind], rr['z'][_iind])
        axtraj2.plot(rr['R'][_iind]*np.cos(rr['phi'][_iind]),rr['R'][_iind]*np.sin(rr['phi'][_iind]))
        #axtraj.plot(rr['R'][_iind][0], rr['z'][_iind][0], 'kx')
        #axtraj.plot([i['R'],rr['R'][_iind][0]] , [i['z'],rr['z'][_iind][0]], 'k-')

    bb=run.bfield.read()
    _Rb = np.linspace(bb['R_min'], bb['R_max'], bb['n_R'])
    _zb = np.linspace(bb['z_min'], bb['z_max'], bb['n_z'])
    #axtraj.contour(_Rb, _zb, bb['psi'].T, bb['psisepx']) 
    axtraj.grid('on')
    axtraj.axis('equal'); axtraj2.axis('equal')
    limit_labels(axtraj, r'R [m]', r'z [m]'); fig_traj.tight_layout()


for ii, run in enumerate(runs):
    e=run.endstate.read()
    i=run.inistate.read()
        
    div=np.nanmax(np.abs(e["B_R"]/e["R"] + e["B_R_dR"] + e["B_phi_dphi"]/e["R"] + e["B_z_dz"]))
    print('divergence: ', div)
    
    ind_wall = np.where(e['endcond']!=10)[0]
    ind_wall = np.where(e['endcond']==8)[0] #wall collision endstate
    e_end = run.endstate['energy']/1.602e-19
    
    print('Lost power to wall in case '+str(ii)+' = '+str(np.sum(e_end[ind_wall]*1.602e-19*e['weight'][ind_wall])*1e-6))
    print('Lost markers =', np.size(ind_wall)/np.size(e['endcond']))
    
    e_ini = run.inistate['energy']/1.602e-19
    print('P inj in case '+str(ii)+' = '+str(np.sum(i['weight']*e_ini*1.602e-19*1e-6)))
    
    v_ini = np.sqrt(i['vR']**2+i['vz']**2+i['vphi']**2)
    pitch = i['vpar']/v_ini
    i['pitch']=pitch
    
    # plot trajectories
    try:
        rr = run.orbits.read()
        plot_trajectory(rr,e,ind_wall, run)
    except AttributeError:
        plot_rz(i,e, ind_wall)
        plot_rhopitch(i, ind_wall)
        plot_phitheta(e, ind_wall)
        plot_phi(e, ind_wall)
        plot_E(run.endstate['energy'], ind_wall)
        
    # time = e['time']
    # ax_etime.scatter((e_end[ind_wall]-e_ini[ind_wall])*1e-3, time[ind_wall], label=label[ii])
    
plt.show()




