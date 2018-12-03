import numpy as np
import matplotlib.pyplot as plt
import ascot_utils
import ascot_particles

eq=input('Which eq?')
fname='/mnt/N41a/home/vallar/ASCOT/runs/runs_lorenzo/FILD_ascot_July18/FILD_ascot/eq{:1d}/ascot_dist.h5'.format(eq)

output_dir ='/mnt/N41a/home/vallar/ASCOT/runs/runs_lorenzo/FILD_ascot_July18/FILD_ascot/results/eq{:1d}/'.format(eq)

print fname
print output_dir

p=ascot_particles.TCV_iniend(fname)
p.plot_histo_wall()
f=plt.gcf()
f.savefig(output_dir+'eq{:1d}_the_wall.png'.format(eq), format='png')
f.savefig(output_dir+'eq{:1d}_the_wall.eps'.format(eq), format='eps')
plt.close()
f=plt.gcf()
f.savefig(output_dir+'eq{:1d}_phi_wall.png'.format(eq), format='png')
f.savefig(output_dir+'eq{:1d}_phi_wall.eps'.format(eq), format='eps')
plt.close()
f=plt.gcf()
f.savefig(output_dir+'eq{:1d}_wall.png'.format(eq), format='png')
f.savefig(output_dir+'eq{:1d}_wall.eps'.format(eq), format='eps')
plt.close()

p.plot_histo_initial_Emax(dE=0.5, max_mark=500)
ans=raw_input('Are markers ok for R? (Y/n)')
if ans=='n':
    max_mark=input('How many markers?')
    p.plot_histo_initial_Emax(dE=0.5, max_mark=max_mark)
f=plt.gcf()
f.savefig(output_dir+'eq{:1d}_R_dE.png'.format(eq), format='png')
f.savefig(output_dir+'eq{:1d}_R_dE.eps'.format(eq), format='eps')
plt.close()
plt.close('all')
p.plot_histo_initial_Emax(max_mark=500)
ans=raw_input('Are markers ok for dE? (Y/n)')
if ans=='n':
    max_mark=input('How many markers?')
    p.plot_histo_initial_Emax(max_mark=max_mark)
f=plt.gcf()
plt.close()
f=plt.gcf()
f.savefig(output_dir+'eq{:1d}_dE.png'.format(eq), format='png')
f.savefig(output_dir+'eq{:1d}_dE.eps'.format(eq), format='eps')
plt.close()



