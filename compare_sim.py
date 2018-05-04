# Script to compare two ascot simulations
# M. Vallar 05/2018

import sys,os
import ascot_distributions
import matplotlib.pyplot as plt
import numpy as np
styles = ['-', '--', '.-']
def main(argv):
    n_timesnaps = (len(argv)-1)/2
    d = np.empty(n_timesnaps, dtype=object)
    fnames = argv[1:n_timesnaps+1]; t=argv[n_timesnaps+1:]
    for i,fname in enumerate(argv[1:n_timesnaps+1]):
        d[i] = ascot_distributions.TCV_1d(fname)
    f=plt.figure(); axcd = f.add_subplot(111)
    f=plt.figure(); axn = f.add_subplot(111)
    f=plt.figure(); axpress = f.add_subplot(111)
    f=plt.figure(); axp = f.add_subplot(111)

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

        d[i].plot_totalpower(axp)
        axp.lines[i].set_label("t="+str(t[i])+" s")
        axp.lines[i].set_linestyle(styles[i])    
    
    axcd.legend(loc='best'); axcd.set_title('Induced shielded current')
    axn.legend(loc='best'); axn.set_title('Density')
    axpress.legend(loc='best'); axpress.set_title('Pressure')
    axp.legend(loc='best'); axp.set_title('Deposited power')
    plt.show()

if __name__=='__main__':
    main(sys.argv)
