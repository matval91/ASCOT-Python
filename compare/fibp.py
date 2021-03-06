"""
Script to compare fibp and evaluate error
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import math

colors=['b','r','g','k','c','m']
col_RGB = np.array([[0, 0, 255], [255, 0, 0], [0, 255, 0], [0,0,0]])
colors_style=['b--','r--','g--','k--','c--','m--']

plt.rc('font', weight='normal')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
plt.rc('figure', facecolor='white')
plt.rc('legend', fontsize=20)

lw=2.3 #line width for plots


# Data should be in the format [rho, fibp0, fibp1,...] with one line of header
data = np.loadtxt('fibp_comp.dat', dtype=float, skiprows=1, unpack=True)
infile=open('fibp_comp.dat')
lines=infile.readlines()
head=lines[0]
labels=head.split()[2:]
labels=[r'$10^{5}$ ptcls',r'$2 \times 10^{5}$ ptcls',r'$5 \times 10^{5}$ ptcls',r'$10^{6}$ ptcls']
infile.close()

#data=np.transpose(data)
rho=data[0,:]
fibp=data[1:,:]
n_fibp=np.shape(fibp)[0]
coeff   = np.zeros((6, n_fibp), dtype=float)
polyval = np.zeros(6, dtype=np.poly1d)
#INTERPOLATION
for i in range(n_fibp):
    tmp=np.polyfit(rho[5:-10], fibp[i][5:-10], 5) 
    coeff[:,i]=tmp
    polyval[i]=np.poly1d(tmp)


RMS=np.zeros(n_fibp,dtype=float)
nptcls = [1e5,2e5,5e5,1e6]
#MEAN ABS DIFFERENCE BTW FIT AND DATA
for i in range(n_fibp):
    tmp=np.mean(np.abs(fibp[i,5:-10]-polyval[i](rho[5:-10])))
    RMS[i]=tmp

# THE ERRORS SHOULD GO AS N^{-1/2}, I fit the squared
RMS2 = 1/np.power(RMS,2)
print RMS
tmp=np.polyfit(np.log(nptcls), np.log(RMS), 1)
error_fit=np.poly1d(tmp)
#x_fun = np.linspace(nptcls[0], nptcls[-1],1000)
#f = lambda x: 1/math.sqrt(x)*RMS[0]/(1/math.sqrt(nptcls[0]))



#PLOTS
f  = plt.figure()
ax = f.add_subplot(111)
for i in range(n_fibp):
    ax.plot(rho[2:-5], fibp[i,2:-5], colors[i], label=labels[i], linewidth=lw)
    ax.plot(rho[5:-10], polyval[i](rho[5:-10]), colors_style[i], linewidth=lw)

ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'Fast Ions Birth profile $\frac{1}{m^{3}s}$')
M = 4
yticks = ticker.MaxNLocator(M)
xticks = ticker.MaxNLocator(M)
# Set the yaxis major locator using your ticker object. You can also choose the minor
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
ax.xaxis.set_major_locator(xticks)
plt.legend(loc='best')

plt.rc('font', weight='normal')
plt.rc('axes', labelsize=15, labelweight='normal', titlesize=24)
a = plt.axes([0.65, 0.275, .2, .2])#, facecolor='y')
#a.plot(nptcls, RMS2,'o')
#a.plot(nptcls, error_fit(nptcls),'k-')
for i in range(len(nptcls)):
    a.plot(np.log(nptcls[i]), np.log(RMS[i]),'o', c=col_RGB[i]/255.)
a.plot(np.log(nptcls), error_fit(np.log(nptcls)),'k-')
#a.plot(np.log(x_fun), np.log(map(f,x_fun)),'k-')

#a.set_xscale('log')
#a.set_yscale('log')
#plt.xlim(0.5e5,1.05e6)
#plt.ylim(0,1e18)
a.set_xlabel('N. markers'); a.set_ylabel(r'$\sigma$ (log scale)')
M = 4
yticks = ticker.MaxNLocator(M)
xticks = ticker.MaxNLocator(M)
# Set the yaxis major locator using your ticker object. You can also choose the minor
# tick positions with set_minor_locator.
a.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
a.xaxis.set_major_locator(xticks)
a.xaxis.grid(which="major")
a.set_yticklabels([])
a.set_xticks(np.log([1e5, 2e5, 5e5, 1e6]))
a.set_xticklabels([r'$10^5$', r'$2\times 10^5$',r'$5\times 10^5$',r'$10^6$'], fontsize=12)
ax.grid('on')
plt.tight_layout()
plt.show()
plt.legend(loc='best')
