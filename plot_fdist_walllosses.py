# Script to plot losses to the wall and 2D spatial Fdist

#M. Vallar - matteo.vallar@igi.cnr.it - 05/2018
import sys,os
import ascot_distributions, ascot_particles
import matplotlib.pyplot as plt

def main(argv):
    fname = argv[1]
    d = ascot_distributions.frzpe(fname)
    p = ascot_particles.h5_particles(fname) 
    d.plot_space()
    ax = plt.gca()
    p.plot_wall_losses(ax)


if __name__=='__main__':
    main(sys.argv)
