# This script is useful to plot the families hitting the wall, coming from PT

import ascot_orbits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import ticker
from ascot_utils import _plot_RZsurf

col = ['k', 'r', 'b', 'g', 'm', 'c']
n_areas=6
scenario=2
if scenario==2:
    fname = '/home/vallar/ASCOT/runs/JT60SA/002/orbits_002037.dat'
    fname_surf = '/home/vallar/ASCOT/runs/JT60SA/002/ascot_002037.h5'
    n_areas=6

    coord_areas = np.array([[[ 1.39228675, -0.8768595 ],
        [ 2.38041434, -1.18264463]],

       [[-1.96135841, -1.80661157],
        [-3.12914557, -1.67024793]],

       [[-1.96135841, -1.80661157],
        [ 0.55387546, -1.85619835]],

       [[ 0.58381872, -1.83140496],
        [ 3.09905259, -1.89752066]],

       [[-3.09920231, -1.84793388],
        [ 0.53890383, -1.91404959]],

       [[-3.08423068, -1.99669421],
        [ 3.09905259, -2.10413223]]])

o = ascot_orbits.SA_orbits(fname, fname_surf = fname_surf)
ind = np.where(o.data_e['endcond']== 3)[0] #wall
r = o.data_e['R'][ind]
z = o.data_e['z'][ind]        
R0 = o.infile['misc/geomCentr_rz'][0]
z0 = o.infile['misc/geomCentr_rz'][1]
theta = np.arctan2(z-z0,r-R0)
phi = o.data_e['phi'][ind]

rho = o.data_i['rho'][ind]
pitch = o.data_i['pitch'][ind]
R = o.data_i['R'][ind]
z = o.data_i['z'][ind]
partdict = o.partdict[ind]

o.plot_histo_wall()
ax=plt.gca()
#plot of rho, phi space
f_rhopitch = plt.figure(figsize=(15,6))
ax_rhopitch = f_rhopitch.add_subplot(121)
ax_rz = f_rhopitch.add_subplot(122)

#plot of orbit
f_orb, ax_orb = plt.subplots(1,n_areas, figsize=(20,7), sharey=True, sharex=True)
for i in range(n_areas):
    _plot_RZsurf(o.Rsurf, o.zsurf, o.RZsurf, ax_orb[i])
    ax_orb[i].plot(o.R_w, o.z_w, lw=2.3, c='k')
    #==============================================
    # SET TICK LOCATION
    #==============================================
    # Create your ticker object with M ticks
    M = 3
    # tick positions with set_minor_locator.
    #ax_orb[i].yaxis.set_major_locator(ticker.MaxNLocator(M))
    #ax.yaxis.set_minor_locator(yticks_m)
    ax_orb[i].xaxis.set_major_locator(ticker.MaxNLocator(M))
    ax_orb[i].set_xlabel('R [m]'); 
    ax_orb[i].axis('equal'); ax_orb[i].grid('on')

ax_orb[0].set_ylabel('Z [m]');


# Now choose the 6 different areas
#coord_areas = np.zeros((n_areas,2,2), dtype=float)
for i in range(n_areas):
    #p0, p1 = plt.ginput(n=2)
    #print(p0, p1)
    #coord_areas[i,0,:] = p0
    #coord_areas[i,1,:] = p1
    p0 = coord_areas[i,0,:] 
    p1 = coord_areas[i,1,:]

    min_x = min(p0[0], p1[0]); max_x = max(p0[0], p1[0])
    min_y = min(p0[1], p1[1]); max_y = max(p0[1], p1[1])
    w = np.abs(p1[0]-p0[0])
    h = np.abs(p1[1]-p0[1])
    r = patches.Rectangle([min_x, min_y], w,h, angle=0.0, fill=0, lw=2.4, color=col[i])
    ax.add_patch(r)
    
    ind = np.where(np.logical_and(
        np.logical_and(phi<max_x, phi>min_x),
        np.logical_and(theta<max_y, theta>min_y)))[0]
    
    ax_rhopitch.scatter(rho[ind], pitch[ind], s=20, color=col[i])
    ax_rz.scatter(R[ind], z[ind], s=20, color=col[i])

    p2plot = partdict[ind[0]]
    ind2plot = p2plot['id'][0]-1

    o._plot_traj_RZ(ax_orb[i], ind2plot, col=col[i])
    
    #o._plot_trajectory(p2plot=p2plot, col=col[i])


#==============================================
# SET TICK LOCATION
#==============================================
# Create your ticker object with M ticks
M = 5
# tick positions with set_minor_locator.
ax_rz.yaxis.set_major_locator(ticker.MaxNLocator(M))
#ax.yaxis.set_minor_locator(yticks_m)
ax_rz.xaxis.set_major_locator(ticker.MaxNLocator(M))
#==============================================
#==============================================
# SET TICK LOCATION
#==============================================
# Create your ticker object with M ticks
M = 5
# tick positions with set_minor_locator.
ax_rhopitch.yaxis.set_major_locator(ticker.MaxNLocator(M))
#ax.yaxis.set_minor_locator(yticks_m)
ax_rhopitch.xaxis.set_major_locator(ticker.MaxNLocator(M))
#==============================================
ax_rhopitch.set_xlabel(r'$\rho$');ax_rhopitch.set_ylabel(r'$\xi (v_\parallel/v)$');
ax_rhopitch.grid('on')
ax_rhopitch.axvline(1., c='k', lw=2.3)
ax_rhopitch.set_ylim([-0.9, -0.3])
ax_rz.set_xlabel(r'R[m]');ax_rz.set_ylabel(r'z [m]');
ax_rz.grid('on')
plt.setp(ax_rz.get_yticklabels()[0], visible=False) 
plt.setp(ax_rhopitch.get_yticklabels()[0], visible=False) 
ax_orb[0].set_xlim([1.5, 4.5])
f_rhopitch.tight_layout()
f_orb.tight_layout()
plt.show()
