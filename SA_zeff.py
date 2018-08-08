"""
Script to create the zeff starting from extreme ideal cases

y2: parabola with maximum at mag axis, minimum (1) at the edge

y8: sixth-grade with minimum on axis (1.5)

"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import ascot_prof


plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
plt.rc('figure', facecolor='white')
plt.rc('legend', fontsize=20)

infile=h5py.File('/home/vallar/ASCOT/runs/JT60SA/005/ascot_005065.h5', 'r+')
vol = infile['distributions/rhoDist/shellVolume'][:]
rho=np.linspace(0,1,len(vol))
f=plt.figure()
ax=f.add_subplot(111)


param2 = [-2.22, 0.0, 3.3]
y2=np.polyval(param2, rho)
mean=np.trapz(y2*vol)/np.trapz(vol)

print(mean)


param_core = [1.6, 0., 0., 0., 0., 0., 1.53]
y8_a=np.polyval(param_core, rho[:-19])
a=np.array([[ 0.88004032,  2.29010417],
       [ 0.89076613,  2.36666667],
       [ 0.90149194,  2.4359375 ],
       [ 0.91362903,  2.5125    ],
       [ 0.92717742,  2.57083333],
       [ 0.94185484,  2.58541667],
       [ 0.9616129 ,  2.49427083],
       [ 0.97346774,  2.37760417],
       [ 0.98419355,  2.26822917],
       [ 0.99604839,  2.17708333]])
param_edge = np.polyfit(a[:,0], a[:,1],8)
#param_edge=np.array([ -2.47658145e+09,   1.85607410e+10,  -6.08418640e+10,
#         1.13935454e+11,  -1.33316296e+11,   9.98104998e+10,
#        -4.66915787e+10,   1.24782240e+10,  -1.45859784e+09])
y8_b = np.polyval(param_edge, rho[-19:])
y8=np.concatenate([y8_a, y8_b])

mean=np.trapz(y8*vol)/np.trapz(vol)
print(mean)

y=np.full(len(vol), 2., dtype=float)
ax.plot(rho,y, lw=2.3, label='Ref.', color='k')
ax.plot(rho, y2, lw=2.3, label='Acc.', color='b')
ax.plot(rho, y8, lw=2.3, label='Edge', color='r')
# Create your ticker object with M ticks
M = 5
yticks = ticker.MaxNLocator(M)
xticks = ticker.MaxNLocator(M)
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
ax.xaxis.set_major_locator(xticks)
#==============================================
ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$Z_{eff}$')
ax.set_ylim([1.2, 4.])
ax.grid('on')
#Removing first point of y-axis
plt.setp(ax.get_yticklabels()[0], visible=False) 
ax.legend(loc='best')
f.tight_layout()


f=plt.figure(); ax=f.add_subplot(111)
dd={'ref':{'z':y,'lab':'Ref.', 'lc':'k'},
'acc':{'z':y2,'lab':'Acc.','lc':'b'},
'edge':{'z':y8,'lab':'Edge', 'lc':'r'}
}
for k in ['ref', 'acc', 'edge']:
    zeff=dd[k]['z']
    p=ascot_prof.SA_datfiles('/home/vallar/JT60-SA/005/input_from_EUrepository/profscenario5.txt', 51, 2, [2,12],[1,6],5,zeff)
    ax.plot(p.rho, p.ni[1,:],  color=dd[k]['lc'], lw=2.3, label=dd[k]['lab'])

ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$n_C [1/m^3]$')
#ax.set_ylim([1.2, 4.])
# Create your ticker object with M ticks
M = 5
yticks = ticker.MaxNLocator(M)
xticks = ticker.MaxNLocator(M)
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
ax.xaxis.set_major_locator(xticks)
#==============================================
ax.grid('on'); ax.set_xlim([0., 1.2])
plt.setp(ax.get_yticklabels()[0], visible=False) 
ax.legend(loc='best')
f.tight_layout()
plt.show()
