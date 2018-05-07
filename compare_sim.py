# Script to compare two ascot simulations
# M. Vallar 05/2018

import sys,os
import ascot_distributions, ascot_prof
import matplotlib.pyplot as plt
import numpy as np
styles = ['-', '--', '.-']

def plot_profiles(p,t):
    p[0].plot_profiles()
    f_prof=plt.gcf()
    p[1].plot_profiles(f=f_prof)
    for ax in f_prof.axes:
        for i, ll in enumerate(ax.lines):
            if len(ax.lines)<=2:
                ll.set_linestyle(styles[i])
                ll.set_label("t="+str(t[i])+" s")
            else:
                if i<2:
                    ll.set_linestyle(styles[0])
                else:
                    ll.set_linestyle(styles[1])
    # Put a legend to the right of the current axis
    f_prof.axes[2].legend(loc='upper center', bbox_to_anchor=(0.1, 1.45))
    f_prof.axes[3].legend_.remove()

def main(argv):
    n_timesnaps = (len(argv)-1)/2
    d = np.empty(n_timesnaps, dtype=object)
    p = np.empty(n_timesnaps, dtype=object)
    fnames = argv[1:n_timesnaps+1]; t=argv[n_timesnaps+1:]
    for i,fname in enumerate(argv[1:n_timesnaps+1]):
        d[i] = ascot_distributions.TCV_1d(fname)
        print('\n Time='+str(t[i])+' s')
        d[i].print_scalars()

        p[i] = ascot_prof.h5_profiles(fname, 51)

    plot_profiles(p,t)

    f=plt.figure(); axcd = f.add_subplot(111)
    f=plt.figure(); axn = f.add_subplot(111)
    f=plt.figure(); axpress = f.add_subplot(111)
    f=plt.figure(figsize=(16,8)); axpe = f.add_subplot(121); 
    axpi=f.add_subplot(122)
    
    for i in range(n_timesnaps):
        d[i].plot_current(axcd)
        axcd.lines[i].set_label("t="+str(t[i])+" s")
        axcd.lines[i].set_linestyle(styles[i])
        
        d[i].plot_n(axn)
        axn.lines[i].set_label("t="+str(t[i])+" s")
        axn.lines[i].set_linestyle(styles[i])

        d[i].plot_p(axpress)
        axpress.lines[i].set_label("t="+str(t[i])+" s")
        axpress.lines[i].set_linestyle(styles[i])

        d[i].plot_pe(axpe)
        axpe.lines[i].set_label("t="+str(t[i])+" s")
        axpe.lines[i].set_linestyle(styles[i])    

        d[i].plot_pi(axpi)

    axpi.lines[0].set_label("i1 t="+str(t[0])+" s")
    axpi.lines[1].set_label("i2 t="+str(t[0])+" s")
    axpi.lines[2].set_label("i1 t="+str(t[1])+" s")
    axpi.lines[3].set_label("i2 t="+str(t[1])+" s")
    
    axpi.lines[0].set_linestyle(styles[0])
    axpi.lines[1].set_linestyle(styles[0])
    axpi.lines[2].set_linestyle(styles[1])
    axpi.lines[3].set_linestyle(styles[1])    
        
    axcd.legend(loc='best'); axcd.set_title('Induced shielded current')
    axn.legend(loc='best'); axn.set_title('Density')
    axpress.legend(loc='best'); axpress.set_title('Pressure')
    axpe.legend(loc='best'); axpe.set_title('Deposited power - electrons'); 
    axpi.legend(loc='best'); axpi.set_title('Deposited power - ions')
    ylim = [0,1800]
    axpe.set_ylim(ylim); axpi.set_ylim(ylim) 

    plt.show()

if __name__=='__main__':
    main(sys.argv)

