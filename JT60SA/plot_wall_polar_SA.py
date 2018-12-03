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

wall = np.loadtxt('/home/vallar/JT60-SA/PARETI_2D_SA/2D_div.txt')
R_w = wall[:,0]
z_w = wall[:,1]
f=plt.figure(figsize=(6,9))
ax=f.add_subplot(111); ax.axis('equal')
ax.plot(R_w, z_w, lw=2.3, c='k')

ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
ax.set_xlim([1.5, 5])
style="Simple,tail_width=1,head_width=8,head_length=10"
kw = dict(arrowstyle=style, color="k")
R0=2.96; d=0.5
x=np.linspace(R0-d,R0+d,1000)
ax.plot(x, np.sqrt(d**2-(x-R0)**2), lw=2.3, c='k')
a = patches.FancyArrowPatch((R0-d, 0.01), (R0-d,0), **kw)
C_big_phi =plt.Circle((4.5, -1.8), 0.2, color='b', fill=False, linestyle='-', linewidth=3.)
x=np.linspace(4.3, 4.7, 100)
ax.plot(x, x-4.5-1.8, ls='-',lw=2.3, c='b')
ax.plot(x, -x+4.5-1.8,  ls='-',lw=2.3, c='b')

ax.add_patch(a)
ax.add_artist(C_big_phi)
ax.text(2.5, 0.5, r'$\theta$', fontsize=48)
ax.text(3.9, 0.1, r'$\theta=0$', fontsize=48)
ax.text(4.4, -1.4, r'$\phi$', fontsize=48, color='b')

ax.axvline(2.96, c='k', ls='--', lw=2.3)
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
