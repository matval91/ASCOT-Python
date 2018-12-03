import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import numpy as np
plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
plt.rc('figure', facecolor='white')
plt.rc('legend', fontsize=20)

wall = np.loadtxt('/home/vallar/TCV/TCV_FW_coord.dat')
R_w = wall[:,0]
z_w = wall[:,1]
f=plt.figure(figsize=(6,9))
ax=f.add_subplot(111); ax.axis('equal')
ax.plot(R_w, z_w, lw=2.3, c='k')

ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
ax.set_xlim([0.45, 1.8])
style="Simple,tail_width=1,head_width=8,head_length=10"
kw = dict(arrowstyle=style, color="k")
R0=0.88; d=0.2
x=np.linspace(R0-d,R0+d,1000)
ax.plot(x, np.sqrt(d**2-(x-R0)**2), lw=2.3, c='k')
a = patches.FancyArrowPatch((R0-d, 0.01), (R0-d,0), **kw)
C_big_phi =plt.Circle((1.4, -0.7), 0.1, color='b', fill=False, linestyle='-', linewidth=3.)
x=np.linspace(1.3, 1.5, 100)
ax.plot(x, x-1.4-0.7, ls='-',lw=2.3, c='b')
ax.plot(x, -x+1.4-0.7,  ls='-',lw=2.3, c='b')

ax.add_patch(a)
ax.add_artist(C_big_phi)
ax.text(0.7, 0.2, r'$\theta$', fontsize=40)
ax.text(1.2, 0.03, r'$\theta=0$', fontsize=40)
ax.text(1.35, -0.5, r'$\phi$', fontsize=40, color='b')

ax.axvline(0.88, c='k', ls='--', lw=2.3)
ax.plot(ax.get_xlim(), [0, 0], c='k', ls='--', lw=2.3)
#ax.grid('on', alpha=0.6)

# Create your ticker object with M ticks
M = 5
yticks = ticker.MaxNLocator(M)
xticks = ticker.MaxNLocator(M)
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
ax.xaxis.set_major_locator(xticks)
#==============================================

f.tight_layout()
plt.show()
