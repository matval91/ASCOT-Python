import ascot_distributions
import ascot_particles
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np

plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
# Create your ticker object with M ticks
M = 4
yticks = ticker.MaxNLocator(M)
yticks_m=ticker.MaxNLocator(M*2)
xticks = ticker.MaxNLocator(M)



styles = ['-','--', '-.']
plot_flag=1
flag_adas=0
mach = 'JT60SA'
mach='TCV'
shot = 2
shot = 53454
ascot_run = 901
if shot==2:
    runs = [22,25,27] #for scenario 2
else:
    runs = [13,15,17]
    runs = [1,5]
shine_thr = [0,0,0]
power = [0,0,0]
en = [500, 350, 250]

for i,e in enumerate(en):
    power[i]=5e6*(e*1000/500000.)**(2.5)
fname = '/home/vallar/ASCOT/runs/TCV/'+"{:05d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(ascot_run)+'.h5'
dis = ascot_distributions.distribution_1d(fname)
len(dis.volumes)
fibp = np.zeros((len(runs), len(dis.volumes)), dtype=float)
for i, el in enumerate(runs):
    dis.TCV_calc_FIBP(0, shot, el)
    fibp[i,:] = dis.fibp
    print "Generated particles: ", dis.fibp_particles
    fname_bb = '/home/vallar/ASCOT/runs/TCV/'+"{:03d}".format(shot)+'/bbnbi_'+"{:03d}".format(shot)+"{:03d}".format(el)+'.h5'

    tmp_ptcls = ascot_particles.h5_particles(fname_bb)
    tmp_ptcls.TCV_calc_shinethr()
    #shine_thr[i] = tmp_ptcls.shinethr_abs[0]


print shine_thr, " W"
for i,el in enumerate(runs):
    print shine_thr[i]/power[i], "%"

if plot_flag == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        labels = str(runs[i])+' run'
        labels = str(en[i])+' keV'
        ax.plot(dis.rho, fibp[i,1:], linewidth = 2.5, linestyle=styles[i],color='k', label = labels)

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'Fast ion birth profile $[1/(m^3 s)]$')
    ax.legend(loc='best')

if plot_flag == 1 & flag_adas==1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        ax.plot(dis.rho, fibp[i,:], linewidth = 3, linestyle=styles[i],color='k', \
                label =  adas_label[i])

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'Fast ion birth profile $[1/(m^3 s)]$')
    ax.legend(loc='best')
        
plt.subplots_adjust(top=0.95,bottom=0.12,left=0.18,right=0.95)
# Set the yaxis major locator using your ticker object. You can also choose the minor
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
ax.xaxis.set_major_locator(xticks)
ax.set_ylim(0,3.0e18)

plt.show()
