import ascot_distributions
import ascot_particles
import matplotlib.pyplot as plt
import numpy as np

styles = ['-','--', '-.']
plot_flag=1

mach = 'TCV'
mach = 'JT60SA'
shot = 29475
shot=002
ascot_run = 12
runs = [3,5]
runs = [30, 32]
labels = ['ADAS','Suzuki']

fname = '/home/vallar/ASCOT/runs/'+mach+'/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(ascot_run)+'.h5'
dis = ascot_distributions.distribution_1d(fname)
fibp = np.zeros((len(runs), len(dis.volumes)), dtype=float)
for i, el in enumerate(runs):
    dis.SA_calc_FIBP(0, shot, el)
    fibp[i,:] = dis.fibp

if plot_flag == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        ax.plot(dis.rho, fibp[i,1:], linewidth = 3, linestyle=styles[i],color='k', label = labels[i])

    ax.legend(loc='best')
    plt.show()
