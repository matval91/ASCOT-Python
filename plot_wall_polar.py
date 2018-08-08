import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

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
f.tight_layout()
plt.show()
