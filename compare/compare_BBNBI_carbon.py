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
shot = 2
shot = 5
ascot_run = 0
if shot == 2:
    runs = [9, 38] #for scenario 2
else:
    runs = [5,20]
shine_thr = [0,0,0]
power = [0,0,0]
carbon= ["No imp.", "Imp."]

fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(ascot_run)+'.h5'
dis = ascot_distributions.distribution_1d(fname)
fibp = np.zeros((len(runs), len(dis.volumes)), dtype=float)
fibp_NNB = np.zeros((len(runs), len(dis.volumes)), dtype=float)
fibp_PNB = np.zeros((len(runs), len(dis.volumes)), dtype=float)

for i, el in enumerate(runs):
    dis.SA_calc_FIBP(0, shot, el)
    fibp[i,:] = dis.fibp
    fibp_NNB[i,:] = dis.fibp_NNB
    fibp_PNB[i,:] = dis.fibp_PNB

    print "Generated particles: ", dis.fibp_particles
    fname_bb = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/bbnbi_'+"{:03d}".format(shot)+"{:03d}".format(el)+'.h5'

    tmp_ptcls = ascot_particles.h5_particles(fname_bb)
    tmp_ptcls.calc_shinethr()
    shine_thr[i] = tmp_ptcls.shinethr_abs[0]
    npart = np.sum(tmp_ptcls.data_i['weight'])
    print " Generated particles from inistate: ", npart

print shine_thr, " W"
if plot_flag == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        labels = str(carbon[i])
        ax.plot(dis.rho, fibp[i,1:], linewidth = 2.5, linestyle=styles[i],color='k', label = labels)
        #labels = str(carbon[i] + ' NNBs')
        #ax.plot(dis.rho, fibp[i,:-1], linewidth = 2.5, linestyle=styles[i],color='r', label = labels)
        #labels = str(carbon[i] + ' PNBs')
        #ax.plot(dis.rho, fibp_PNB[i,:-1], linewidth = 2.5, linestyle=styles[i],color='b', label = labels)
        
    ax.set_ylim([0,2.2e19])
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'Fast ion birth profile $[1/(m^3 s)]$')
    ax.legend(loc='best')
    
    plt.subplots_adjust(top=0.95,bottom=0.12,left=0.18,right=0.95)
    # Set the yaxis major locator using your ticker object. You can also choose the minor
    # tick positions with set_minor_locator.
    ax.yaxis.set_major_locator(yticks)
    #ax.yaxis.set_minor_locator(yticks_m)
    ax.xaxis.set_major_locator(xticks)
    ax.set_ylim(0,4e19)

    plt.show()
