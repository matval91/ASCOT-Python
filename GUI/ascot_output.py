import numpy as np
#GUI
import Tkinter as tk
# for the definition of appropriate plotting style
import tksty, tkFileDialog
import h5py
import matplotlib.pyplot as plt
from matplotlib import ticker
import os 

import ascot_particles
import ascot_prof
import ascot_Bfield
import ascot_distributions

sty = tksty.TKSTY()

def main():
    global ascot_fname
    global dis
    global prof
    global bfield
    global part_b

    
    #tmp_fname = tkFileDialog.Open(filetypes=[('hdf5', '*.h5'),('All',''*)],\
    #                              initialdir=os.environ['HOME'], title='Open File ASCOT').show()
    tmp_fname = tkFileDialog.askopenfilenames(title='Open ASCOT and BBNBI files')
    if 'ascot_' in tmp_fname[0]:
        ind_asc=0
        ind_bb=1
    else:
        ind_asc=1
        ind_bb=0
    ascot_fname = tmp_fname[ind_asc]

    if tmp_fname:    
        dis    = ascot_distributions.distribution(tmp_fname[ind_asc])
        prof   = ascot_prof.h5_profiles(tmp_fname[ind_asc])
        bfield = ascot_Bfield.Bfield_out(tmp_fname[ind_asc])
        part_b = ascot_particles.particles(tmp_fname[ind_bb])
        
        if __name__ == '__main__':
            fmframe = tk.TK()
        else:
            fmframe = tk.Toplevel()
    
    
        commframe_f(fmframe)
        bbnbiframe_f(fmframe)
        ascotframe_f(fmframe)
    
        qframe = sty.mySubframe(fmframe)
        sty.myButton(qframe, 'QUIT', fmframe.destroy, bgc=tksty.qtcol)



def commframe_f(fmframe):

    commframe = sty.mySubframe(fmframe)
    #sty.myButton(commframe, 'PROFILES', prof.checkplot  , nrow=1, ncol=1)
    sty.myButton(commframe, 'MAGNETIC', bfield.checkplot, nrow=1, ncol=2)

def bbnbiframe_f(fmframe):
    bbnbiframe = sty.mySubframe(fmframe)
    global store_data
    store_data=tk.IntVar(); store_data.set(0)

    sty.myButton(bbnbiframe, 'Ionisation in XY space', part_b.plot_XYpart  , nrow=1, ncol=1)
    sty.myButton(bbnbiframe, 'Shine Through'      , part_b.calc_shinethr, nrow=1, ncol=2)
    sty.myButton(bbnbiframe, 'Fast Ion birth profile', calc_FIBP, nrow=2, ncol=1)
    sty.myCb(bbnbiframe,txt='Save data?', var=store_data, nrow=2, ncol=2)


def ascotframe_f(fmframe):
    ascotframe = sty.mySubframe(fmframe)
    sty.myButton(ascotframe, 'All distributions', plot_distributions_GUI      , nrow=1, ncol=1)
    sty.myButton(ascotframe, 'Dis w ions'  , plot_distributions_wions, nrow=1, ncol=2)
    sty.myButton(ascotframe, 'Scalar quantities', dis.print_scalars, nrow=2,ncol=1)
    
def plot_distributions_GUI():
    global plot_n, plot_eden, plot_jpar, plot_pel, plot_pi1, plot_torque, plot_ctorque
    plot_n=tk.IntVar(); plot_n.set(1)
    plot_eden=tk.IntVar(); plot_eden.set(1)
    plot_jpar=tk.IntVar(); plot_jpar.set(1)
    plot_pel=tk.IntVar(); plot_pel.set(1)
    plot_pi1=tk.IntVar(); plot_pi1.set(1)
    plot_torque=tk.IntVar(); plot_torque.set(1)
    plot_ctorque=tk.IntVar(); plot_ctorque.set(1) 
    
    if __name__ == '__main__':
        fmframe = tk.TK()
    else:
        fmframe = tk.Toplevel()  
        """
        self.name_dict = {'n':1, 'e_den':2,\
                          'jpar':3, 'jperp':4, \
                          'jxB':5, 'jxBstate':6, \
                          'CX_ionsource':7, 'CX_ionensource':8,\
                          'CX_neutsource':9, 'CX_neutensource':10,\
                          'torque':11,'par_e_den':12,\
                          'tot_tor_j':13, 'ptot':14, 'ppar':15, 'pperp':16,\
                          'flr_torque':17,\
                          'th_n':18, 'th_e_n':19, 'th_torque':20, 'abs_ICRH':21,\
                          'pel':22, 'pi1':23,\
                          'ctor_el':24, 'ctor_i1':25\
                          }
        """
    #density,pressure, jpar, pel, pi1, torque, e_density 
    sty.myCb(fmframe,txt='Density',            var=plot_n, nrow=1, ncol=1)
    sty.myCb(fmframe,txt='Energy den',     var=plot_eden, nrow=1, ncol=2)
    sty.myCb(fmframe,txt='Parallel j',         var=plot_jpar, nrow=2, ncol=1)
    sty.myCb(fmframe,txt='Power to e', var=plot_pel, nrow=2, ncol=2)
    sty.myCb(fmframe,txt='Power to i',      var=plot_pi1, nrow=3, ncol=1)
    sty.myCb(fmframe,txt='Torque (jxB)',       var=plot_torque, nrow=3, ncol=2)
    sty.myCb(fmframe,txt='Coll torque', var=plot_ctorque, nrow=4, ncol=1)
    sty.myButton(fmframe, 'Plot', plot_distributions, nrow=5, ncol=1)
    sty.myButton(fmframe, 'Dismiss', fmframe.destroy, bgc=tksty.qtcol, nrow=5, ncol=2)
 
    
      

def plot_distributions():
    plot_array = []    
    if plot_n      .get()==1: plot_array.append(0)
    if plot_eden   .get()==1: plot_array.append(1)
    if plot_jpar   .get()==1: plot_array.append(2)
    if plot_pel    .get()==1: plot_array.append(21)
    if plot_pi1    .get()==1: plot_array.append(22)
    if plot_torque .get()==1: plot_array.append(4)
    if plot_ctorque.get()==1: plot_array.append(24)  
    
    if len(plot_array)!=0:
        dis.SA_plot_profs(plot_array)
    else:
        dis.SA_plot_profs()
       
    #dis.plot_profs(['n', 'e_den'])
    #dis.plot_profs(['ptot', 'ppar', 'pperp'])
    #dis.plot_profs(['ptot', 'ppar', 'pperp'])

def plot_distributions_wions():
    dis.plot_profs_wions(['pel','pi1','pi2','pi3'])


def calc_FIBP():
    try:
        volumes = dis.volumes
    except:
        dis._evaluate_shellVol()
        volumes = dis.volumes
    global ascot_fname
    infile=h5py.File(ascot_fname)
    rho=infile['distributions/rhoDist/abscissae/dim1'].value
    fibp=np.zeros(len(rho),dtype=float)
    fibp_PNB=np.zeros(len(rho),dtype=float)
    fibp_P_tanj=np.zeros(len(rho),dtype=float)
    fibp_P_perp=np.zeros(len(rho),dtype=float)
    fibp_NNB=np.zeros(len(rho),dtype=float)
    rho_edg = rho+(rho[-1]-rho[-2])*0.5

    origins=part_b.data_i['origin']
    part_rho = part_b.data_i['rho']
    for i,r in enumerate(rho):
        if r==rho[-1]:
            continue
        ind = [(part_rho>rho_edg[i]) &  (part_rho<rho_edg[i+1])]
        ind_P = [(part_rho>rho_edg[i]) & (part_rho<rho_edg[i+1]) & (origins != 3031) & (origins != 3032)]
        ind_P_perp = [(part_rho>rho_edg[i]) & (part_rho<rho_edg[i+1]) & \
                      (origins != 3031) & (origins != 3032) & \
                      (origins != 309) & (origins != 310) & \
                      (origins != 312) & (origins != 311) & \
                      (origins != 3637) & (origins != 3638) & \
                      (origins != 3639) & (origins != 3640)\
                      ]
        
        weight = np.sum(part_b.data_i['weight'][ind])
        weight_P = np.sum(part_b.data_i['weight'][ind_P])
        weight_PP = np.sum(part_b.data_i['weight'][ind_P_perp])

        fibp[i] = weight/volumes[i]
        fibp_PNB[i]=weight_P/volumes[i]
        fibp_P_perp[i]=weight_PP/volumes[i]
        fibp_P_tanj[i]=fibp_PNB[i]-fibp_P_perp[i]
        fibp_NNB[i]=fibp[i]-fibp_PNB[i]

    if store_data.get()==1:
        header='rho,   total, P-perp,   P-tanj,   NNB'
        np.savetxt('fibp.dat', np.transpose([rho,fibp,fibp_P_perp, fibp_P_tanj, fibp_NNB]), header=header, fmt='%1.4e')
        print 'FIBP stored in fibp.dat'
                

##     plt.plot(rho,fibp)
##     plt.xlabel('rho')
##     plt.ylabel('Fast Ions Birth Profile (1/(s*m^3))')
##     plt.show()

    plot_article(4,[rho, fibp_P_tanj,  fibp_P_perp, fibp_NNB, fibp], ['P-T', 'P-P', 'NNB', 'TOT'],r'$\rho$','Fast Ions Birth Profile (1/($s\cdot m^3$))')


def plot_article(n_lines, data, data_labels, xlabel, ylabel):
        #=====================================================================================
        # SET TEXT FONT AND SIZE
        #=====================================================================================
        plt.rc('font', family='serif', serif='Palatino')
#        plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=====================================================================================
        col=['r','g','b','k']

        fig=plt.figure()
        ax=fig.add_subplot(111)
        for i in range(n_lines):
            ax.plot(data[0], data[i+1], label=str(data_labels[i]),color=col[i], linewidth=3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_ylim([0,2.5e19])
        #=====================================================================================
        # ADJUST SUBPLOT IN FRAME
        #=====================================================================================
        plt.subplots_adjust(top=0.95,bottom=0.12,left=0.12,right=0.95)
        #=====================================================================================

        #=====================================================================================
        # SET TICK LOCATION
        #=====================================================================================

        # Create your ticker object with M ticks
        M = 4
        yticks = ticker.MaxNLocator(M)
        xticks = ticker.MaxNLocator(M)
        # Set the yaxis major locator using your ticker object. You can also choose the minor
        # tick positions with set_minor_locator.
        ax.yaxis.set_major_locator(yticks)
        #ax.yaxis.set_minor_locator(yticks_m)
        ax.xaxis.set_major_locator(xticks)
        #=====================================================================================
        ax.legend(loc='best')



        #ax.set_ylim([0,2.5e19])
