import a5py.ascot5io.ascot5 as a5
from utils.plot_utils import _plot_2d, limit_labels
import matplotlib.pyplot as plt
import numpy as np


label = ['2D', 'Ripple']
col = ['k', 'r']
f=a5.Ascot('/home/vallar/ASCOT/runs/SA_003/3D/NNB_3D_ripple/ascot.h5')
# B_3DS-3439715758 has cocos 7. B0 is fixed to Raxis
# B_3DS-3519459860 is fixed so that B0 is at R0. WRONG GRID
# B_3DS-3241863959 should be ok
# B_3DS-2300597066 with no -1 factor in producing 3D bfield
# B_3DS-3632670271
#2D field used is B_2DS-0461682159

w=f.wall.wall_2D_3851912840.read()
Rw=w['R']; zw=w['z']

f1 = plt.figure(); ax1 = f1.add_subplot(111)
f2 = plt.figure(); ax2 = f2.add_subplot(111)
for ii, b in enumerate([f.bfield.B_2DS_0461682159, f.bfield.B_3DS_3632670271]):
    b=b.read()
    bphi = b['B_phi']; br = b['B_R']
    R = np.linspace(b['R_min'], b['R_max'], b['n_R'][0])
    z = np.linspace(b['z_min'], b['z_max'], b['n_z'][0])
    psi = b['psi']
    if len(np.shape(bphi))>2:
        bphi = bphi[15,:,:]
        br = br[0,:,:]
    #CS=ax1.contour(R,z,bphi, [-4, -3., -2., -1.]); plt.colorbar(CS)
    ax1.plot(R, bphi[int(b['n_z'][0]/2),:], label=label[ii])
    try:
        CS=ax2.contour(psi, [-1,0,2,4]); plt.colorbar(CS)
    except:
        continue

ax2.plot(Rw, zw, 'k', lw=2.3)
ax1.legend(loc='best')
plt.show()
    
