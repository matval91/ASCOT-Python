import ascot_particles
import matplotlib.pyplot as plt
import numpy as np

plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)


styles = ['-','--', '-.']
color=['r','b']
shot = 2
runs = [30,31]

npart = 100

f = plt.figure()
ax = f.add_subplot(111)
for i, el in enumerate(runs):
    fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(el)+'.h5'
    prt = ascot_particles.h5_particles(fname)
    x = prt.data_e['time']
    y = prt.data_e['energy']
    ax.plot(x,y, color[i]+'o', ms=10, label="RUN "+str(el))

ax.legend(loc='best')
ax.set_xlabel('Time')
ax.set_ylabel('Energy')
plt.show()
