from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, patches
import collections
import math, platform
import h5py
from utils.plot_utils import _plot_2d, _plot_1d

colours = ['k', 'g','b','r','c']
styles = ['-','--','-.']

cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

my_cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)

class particles:
    """
    SUPERCLASS
    
    """
    def __init__(self, fname_surf):
        # Initialisation
        self.npart   = 0
        self.nfields = 0
        self.field   = []
        self.unit    = []
        self.data_i  = collections.defaultdict(list)
        self.data_e  = collections.defaultdict(list)
        self.R_w = []
        self.z_w = []
        self.Rsurf = []
        self.zsurf = []
        self.RZsurf = []
        self.Rsurf, self.zsurf, self.RZsurf = _RZsurf(fname_surf)
        self.R_w, self.z_w = _readwall(self.device)
        try:
            self.R0 = h5py.File(fname_surf)['misc/geomCentr_rz'][0]
            self.z0 = h5py.File(fname_surf)['misc/geomCentr_rz'][1]
        except:
            self.R0=0
            self.z0=0
            print("Impossible to read R0 and Z0")

    def _calc_weight(self):
        """
        hidden method to calculate the weight of the particles for each origin: in the inistate,
        in the final state and the total one.
        """
        try:
            self.data_i['weight'][:].mean()
        except:
            print("No weights available")
            return
        
        #self._calc_originbeam()
        self.w_i = np.zeros((len(self.origins)), dtype=float)
        self.w_e = np.zeros((len(self.origins)), dtype=float)
        self.w   = np.zeros((len(self.origins)), dtype=float)

        for ind_i, i in enumerate(self.origins):
            ind = self.data_i['origin'][:]==i
            self.w_i[ind_i] = sum(self.data_i['weight'][ind])
            ind = self.data_e['origin'][:]==i
            self.w_e[ind_i] = sum(self.data_e['weight'][ind])  
            self.w[ind_i]  = self.w_i[ind_i]+self.w_e[ind_i]

    def plot_histo_wall(self, ax=0, lastpoint=1):
        """
        Histogram of deposition to the wall
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        R0 = self.infile['misc/geomCentr_rz'][0]
        try:
            z0 = self.infile['misc/geomCentr_rz'][1]
        except:
            z0=0
            print("z0 set to 0")

        theta = np.arctan2(z-z0,r-R0)
        phi = self.data_e['phi'][ind]
        #wallrz= [self.R_w, self.z_w]
        #energy = self.data_e['energy'][ind]*1e-3
        angles_ticks=[-3.14, -1.57, 0., 1.57, 3.14]
        angles_labels=[r'-$\pi$',r'-$\pi/2$',r'0',r'$\pi/2$',r'$\pi$']
        _plot_2d(phi, theta, xlabel=r'$\phi$ [rad]',ylabel=r'$\theta$ [rad]', hist=1, xlim=[-3.14, 3.14],ylim=[-3.14, 3.14], cblabel='# markers', lastpoint=lastpoint)
        ax=plt.gca()
        ax.set_xticks(angles_ticks); ax.set_xticklabels(angles_labels)
        ax.set_yticks(angles_ticks); ax.set_yticklabels(angles_labels)

        _plot_1d(theta, xlabel=r'$\theta$ [rad]', hist=1)
        ax=plt.gca()
        ax.set_xticks(angles_ticks); ax.set_xticklabels(angles_labels)
        _plot_1d(phi, xlabel=r'$\phi$ [rad]', hist=1)
        ax=plt.gca()
        ax.set_xticks(angles_ticks); ax.set_xticklabels(angles_labels)
        ax=plt.gca()
        r = self.data_i['R'][ind]
        z = self.data_i['z'][ind]
        theta = np.arctan2(z-z0,r-R0)
        phi = self.data_i['phi'][ind]
        #_plot_2d(phi, theta, xlabel='', ylabel='', ax=ax, scatter=1)
                                    
class h5_particles(particles):
    """
    superClass (inherited from particles) handling h5 files (e.g. bbnbi.h5, ascot.h5)
    """
    def __init__(self, fname, fname_surf=''):
        """
        Initialising
        """
        self._fname = fname
        indd = self._fname[-9:-3]
        self.id = indd
        if indd[0:2] != '00':
            self.id = self._fname[-11:-3]
        if fname_surf=='':
            fname_surf=fname
        particles.__init__(self, fname_surf)
        dataf = h5py.File(self._fname)
        self.infile=dataf
        if 'shinethr' not in self.infile.keys():
            print(list(dataf.keys()))
        
        try:
            self._endgroup = 'shinethr' #this is for bbnbi
            self.field = dataf[self._endgroup].keys()
            self.origins = np.sort(np.array(list(set(self.data_i['origin']))))
        except:
            self._endgroup = 'endstate' #this is for ascot
            self.field = dataf[self._endgroup].keys()
            self.origins = np.sort(np.array(list(set(self.data_i['origin']))))
                    
        self.nfields=len(self.field)      
        for key in self.field:
            tt=float
            if key=='id':
                tt=int
            self.data_i[key] = np.array(dataf['inistate/'+key], dtype=tt)
            self.data_e[key] = np.array(dataf[self._endgroup+'/'+key], dtype=tt)
                
        # evaluate number of particles
        self.npart = np.max(self.data_i['id'])
            
        #CONVERT PHI FROM DEG TO RAD
        self.data_i['phiprt']=self.data_i['phiprt']*math.pi/180.
        self.data_e['phiprt']=self.data_e['phiprt']*math.pi/180.
        self.data_i['phi']=self.data_i['phi']*math.pi/180.

        #print("FIELDS IN FILE ",self._fname,", section inistate and ",self._endgroup," :")
        #print(list(self.field))
        
        # Swapping R and Rprt, since they are opposite to matlab
        self.data_e['R']   = self.data_e['Rprt']
        self.data_e['z']   = self.data_e['zprt']
        self.data_e['phi'] = self.data_e['phiprt']


    def plot_RZ(self, ax=0, shpart=0):
        """
        Method to plot R vs z of the ionised particles, useful mostly 
        with bbnbi
        """
        try:
            x=self.data_i['Rprt']
            y=self.data_i['zprt']
        except:
            x=self.data_i['R']
            y=self.data_i['z']

        if np.mean(x)==999. or np.mean(x)==-999.:
            x=self.data_i['R']
            y=self.data_i['z']
        xlab = 'R [m]'
        ylab = 'z [m]'
        wallrz= [self.R_w, self.z_w]
        surf=[self.Rsurf, self.zsurf, self.RZsurf]
        _plot_2d(x, y, xlabel=xlab, ylabel=ylab,  title='RZ ionization',\
                 wallrz=wallrz, surf=surf, ax=ax, xlim=self.xlim, scatter=1)


        if 'bbnbi' in self._fname and shpart!=0:
            ax = plt.gca()
            R = self.data_e['R']
            ind = np.where(R<4.5)
            R=R[ind]
            phi = self.data_e['phi'][ind]
            x,y = R*np.cos(phi), R*np.sin(phi)
            z = self.data_e['z'][ind]
            _plot_2d(R,z, scatter=1, ax=ax)

    def plot_XY(self, ax=0, shpart=0):
        """
        Method to plot XY of ionisation, without difference between the beams
        """
        try:
            R=self.data_i['Rprt']
            z=self.data_i['zprt']
            phi = self.data_i['phiprt']
        except:
            R=self.data_i['R']
            z=self.data_i['z']
            phi = self.data_i['phi']

        if np.mean(R)==999. or np.mean(R)==-999.:
            R=self.data_i['R']
            z=self.data_i['z']
            phi = self.data_i['phi']

        x=np.zeros(self.npart)
        y=np.zeros(self.npart)
        x = R*np.cos(phi)
        y = R*np.sin(phi)
        xlab = 'X [m]'
        ylab = 'Y [m]' 
        R0=self.R0
            
        wallxy= [self.R_w, self.z_w]
        
        _plot_2d(x, y, xlabel=xlab, ylabel=ylab,  title='XY Ionization',\
                 wallxy=wallxy, R0=R0, ax=ax)
        if 'bbnbi' in self._fname and shpart!=0:
            ax = plt.gca()
            R = self.data_e['R']
            ind = np.where(R<4.5)
            R=R[ind]
            phi = self.data_e['phi'][ind]
            x,y = R*np.cos(phi), R*np.sin(phi)
            _plot_2d(x,y, scatter=1, ax=ax)

    def endcondition(self):
        """
        Computes the endcondition of the particles and stores a dict with the amount of MC particles with that final condition
        """
        self.ptot_inj  = np.dot(self.data_i['weight'],self.data_i['energy'])*1.602e-19
        self.plost_w   = 0.

        pieplot_label = ['CPU', 'tmax', 'emin', 'wall', 'th.']
        pieplot_x = dict.fromkeys(pieplot_label, 0.)
        npart = float(len(self.data_e['endcond']))

        errendcond  = {'aborted':-2, 'rejected':-1}
        counter = 0; pwall=0; pth=0;
        for key in errendcond:
            ind = np.where(self.data_e['endcond']==errendcond[key])[0]
            if len(ind)!=0:
                counter += len(ind)
                pieplot_x['CPU'] += len(ind)/npart
                print("STRANGE END CONDITION! {} ({:d}) : {:d} ({:.2f}%) particles".format(key, errendcond[key], len(ind), len(ind)/npart*100.))

        physcond = {'none':0,'tmax':1,'emin':2,'wall':3,'th.':4}
        for key in physcond:
            ind = np.where(self.data_e['endcond']== physcond[key])[0]
            counter += len(ind)
            if key!='none':
                pieplot_x[key] += len(ind)/npart
            pow_lost = np.dot(self.data_e['energy'][ind], self.data_e['weight'][ind])*1.602e-19*1e-3
            print("{} ({:2d}) : {:4d} ({:.2f}%) particles, {:7.2f} kW".format(key, physcond[key], len(ind), len(ind)/npart*100., pow_lost))
            if key == 'wall':
                pwall = len(ind)/npart*100.
                self.plost_w = pow_lost*1e3
            elif key == 'th.':
                pth = len(ind)/npart*100.

        infcond = {'cputmax':17, 'outsidemaxrho':10, 'outsidewall':16}
        for key in infcond:
            ind = np.where(self.data_e['endcond']== infcond[key])[0]
            counter += len(ind)
            pieplot_x['CPU'] += len(ind)/npart
            print("{} ({:2d}) : {:4d} ({:.2f}%) particles".format(key, infcond[key], len(ind), len(ind)/npart*100.))    
         
        print("Total particles counted :", counter/npart*100., ' %')
        
        x = np.array(pieplot_x.values())
        #_plot_pie(pieplot_x.values(), lab=pieplot_x.keys(), Id=self.id)

    def plot_initial_wall_pos(self, theta_f=[-3.14, 3.14], phi_f=[-3.14, 3.14]):
        """
        Plot with final position of the particles colliding with wall and energy 
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        xlim = [0, 25.]; ylim=[0., 1.]
        #Filtering on the range desired in phi, theta
        if np.mean(theta_f)==0. and np.mean(phi_f)==0.:
            self.plot_histo_wall()
            ax=plt.gca()
            p1, p2 = plt.ginput(n=2)
            theta_f = [p1[1], p2[1]]
            phi_f = [p1[0], p2[0]]   
            plt.close()
        if theta_f[0]==theta_f[1] and phi_f[0]==phi_f[1]:
            theta_f = [-3.14, 3.14]
            phi_f = [-3.14, 3.14]    
          
        #==========================================================
        # Finding index where particle meet the conditions
        #==========================================================        
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]        
        R0 = self.infile['misc/geomCentr_rz'][0]
        z0 = self.z0
        theta = np.arctan2(z-z0,r-R0)
        phi = self.data_e['phiprt'][ind]
        ind_inist = np.where(np.logical_and(\
                                np.logical_and(\
                                               theta>np.min(theta_f), theta<np.max(theta_f)\
                                           ),\
                                np.logical_and(\
                                               phi>np.min(phi_f), phi<np.max(phi_f)\
                                           )))
        ind_inistate = ind[ind_inist]
        n_wall = len(ind_inistate)
        
        #==========================================================
        # Define inistate variables
        #==========================================================
        r = self.data_i['R'][ind_inistate]
        z = self.data_i['z'][ind_inistate]
        rho = self.data_i['rho'][ind_inistate]
        energy = self.data_i['energy'][ind_inistate]*1e-3
        pitch = self.data_i['pitch'][ind_inistate]
        wallrz= [self.R_w, self.z_w]
        phi = self.data_i['phi'][ind_inistate]
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        ind_lt_r0 = r<R0
        self.plot_iniendstate(theta_f, phi_f, x,y,energy,wallrz, pitch,rho, ind_lt_r0,r,z)

        #==========================================================
        # Define enstate variables
        #==========================================================
        r = self.data_e['R'][ind_inistate]
        z = self.data_e['z'][ind_inistate]
        rho = self.data_e['rho'][ind_inistate]
        energy = self.data_e['energy'][ind_inistate]*1e-3
        deltaE = (self.data_i['energy'][ind_inistate]-\
                  self.data_e['energy'][ind_inistate])*1e-3
        pitch = self.data_e['pitch'][ind_inistate]
        wallrz= [self.R_w, self.z_w]
        phi = self.data_e['phi'][ind_inistate]
        
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        
        self.plot_iniendstate(theta_f, phi_f, x,y,energy,wallrz, pitch,rho, ind_lt_r0,r,z)

        # particles born with R<R0:
        r = self.data_i['R'][ind_inistate]
        ind_r0 = r<self.R0
        r = self.data_i['R'][ind_inistate][ind_r0]
        plost_rltr0 = np.dot(self.data_e['energy'][ind_inistate][ind_r0],self.data_e['weight'][ind_inistate][ind_r0]) 
        self.plost_rltr0 = plost_rltr0*1.602e-19
        try:
            self.plost_w
        except:
            self.endcondition()
        self.plost_rgtr0 = self.plost_w-self.plost_rltr0
        print('PARTICLES WITH R<R0: {:d} = {:.2f} %'.format(len(r), float(len(r))/n_wall*100.))
        # z = self.data_i['z'][ind_inistate][ind_r0]
        # rho = self.data_i['rho'][ind_inistate][ind_r0]
        # energy = self.data_i['energy'][ind_inistate][ind_r0]*1e-3
        # pitch = self.data_i['pitch'][ind_inistate][ind_r0]
        # wallrz= [self.R_w, self.z_w]
        # phi = self.data_i['phi'][ind_inistate][ind_r0]
        deltaE = (self.data_i['energy'][ind_inistate]-\
                  self.data_e['energy'][ind_inistate])*1e-3
        # # Plot of DE of particles born at R<R0
        # fig = plt.figure(figsize=(10,12), dpi=70)
        # fig.canvas.set_window_title('End State '+self.id)
        # fig.text(0.1, 0.01, self.id)
        # axdE = fig.add_subplot(111)
        # _plot_1d([deltaE[~ind_lt_r0], deltaE[ind_lt_r0]], xlabel=r'$\Delta E$ [keV]', ylabel=r'# markers', Id=self.id, hist=1, ax=axdE, multip=1)    


    def plot_iniendstate(self, theta_f, phi_f, x,y,energy,wallrz, pitch,rho, ind_lt_r0,r,z): 
        #==========================================================
        # PLOT OF STATE
        #==========================================================
        fig = plt.figure(figsize=(10,12), dpi=70)
        fig.text(0.1, 0.01, self.id)
        fig.canvas.set_window_title('Initial State '+self.id)
        tit = r'$\theta$=[{:.2f}$\div${:.2f}] | $\phi$=[{:.2f}$\div${:.2f}]'.\
              format(np.min(theta_f), np.max(theta_f), np.min(phi_f), np.max(phi_f))
        fig.suptitle(tit, fontsize=16)
        nrow=3; ncol=2
        #axrz = fig.add_subplot(nrow, ncol, 1)
        #_plot_2d(r, z, 'R [m]', 'z [m]', scatter=energy, Id=self.id, wallrz=wallrz, \
        #         surf=[self.Rsurf, self.zsurf, self.RZsurf], ax=axrz, multip=1)
        #axrz.axvline(x=np.max(self.R_w-0.1))
        #axrz.axvline(x=np.min(self.R_w+0.1))
        axxy = fig.add_subplot(nrow, ncol,2)
        _plot_2d(x, y, 'X [m]', 'Y [m]', scatter=energy, wallxy=wallrz, ax=axxy, multip=1, R0=self.R0)
        axep = fig.add_subplot(nrow,ncol,3)
        _plot_2d(energy, pitch, xlabel=r'E [keV]',ylabel=r'$\xi$', Id = self.id, hist=1, ax=axep, multip=1,\
                 xlim=[0, 30], ylim=[-1.,1.])
        
        axrho = fig.add_subplot(nrow, ncol,4)
        _plot_1d([rho[~ind_lt_r0], rho[ind_lt_r0]], xlabel=r'$\rho$', ylabel=r'# markers', Id=self.id, hist=1, ax=axrho, multip=1)
        axr = fig.add_subplot(nrow, ncol,5)
        _plot_1d([r[~ind_lt_r0], r[ind_lt_r0]], xlabel=r'R [m]', ylabel=r'# markers', Id=self.id, hist=1, ax=axr, multip=1)

        axr.axvline(x=self.R0)
        axr2 = axr.twiny()
        axr2.set_xlim(axr.get_xlim())
        axr2.set_xticks([self.R0])
        axr2.set_xticklabels(['R0']); plt.setp(axr2.get_xticklabels(), rotation='45', fontsize=16)
        axz = fig.add_subplot(nrow, ncol,6)
        _plot_1d([z[~ind_lt_r0], z[ind_lt_r0]], xlabel=r'z [m]', ylabel=r'# markers', Id=self.id, hist=1, ax=axz, multip=1)
        axz.axvline(x=self.z0)
        axz2 = axz.twiny()
        axz2.set_xlim(axz.get_xlim())
        axz2.set_xticks([self.z0])
        axz2.set_xticklabels(['z0']); plt.setp(axz2.get_xticklabels(), rotation='45', fontsize=16)
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)


    def plot_wall_losses(self, ax=0):
        """
        Plot with final position of the particles colliding with wall and energy 
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        R0=self.R0
        wallrz= [self.R_w, self.z_w]
        energy = self.data_e['energy'][ind]*1e-3

        _plot_2d(r,z, 'R [m]', 'z [m]', scatter=energy, wallrz=wallrz, surf=[self.Rsurf, self.zsurf, self.RZsurf], ax=ax)
        

    def plot_histo_rho(self):
        """Plot initial rho of lost particles
        """
        try:
            self.R_w.mean()
        except:
            self._readwall() 
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        x = self.data_i['rho'][ind]
        _plot_1d(x, xlabel=r'$\rho$', ylabel=r'# markers', hist=1, ylim=[0, 2500])
        ax=plt.gca()
        ax.axvline(1.,color='k', lw=2.3)

    def plot_histo_R(self, max_mark=500):
        """Plot initial R of lost particles
        """
        try:
            self.R_w.mean()
        except:
            self._readwall() 

        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        x = self.data_i['R'][ind]
        _plot_1d(x, xlabel=r'R [m]', ylabel=r'# markers', hist=1, ylim=[0, max_mark], color='b')
        ax=plt.gca()
        ax.axvline(self.R0,color='k', lw=2.3)
        rect = patches.Rectangle((max(self.R_w),0),0.05,ax.get_ylim()[-1],edgecolor='k',facecolor='k')
        ax.add_patch(rect)
        rect = patches.Rectangle((min(self.R_w)-0.05,0),0.05,ax.get_ylim()[-1],edgecolor='k',facecolor='k')
        ax.add_patch(rect)
        ax.set_xlim([min(self.R_w)-0.05, max(self.R_w)+0.05])


    def plot_histo_z(self):
        """Plot initial z of lost particles
        """
        try:
            self.R_w.mean()
        except:
            self._readwall() 

        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        x = self.data_i['z'][ind]
        _plot_1d(x, xlabel=r'z [m]', ylabel=r'# markers', hist=1, ylim=[0, 1500])
        ax=plt.gca()
        ax.axvline(self.z0,color='k', lw=2.3)

    def plot_histo_E(self, max_mark=500, ax=0):
        """
        Plot initial E of lost particles
        """
        try:
            self.R_w.mean()
        except:
            self._readwall() 

        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        x_final = self.data_e['energy'][ind]*1e-3
        x_initial = self.data_i['energy'][ind]*1e-3
        x = x_initial-x_final
        x = x[np.where(x<100)[0]]
        _plot_1d(x, xlabel=r'$\Delta E$ [keV]', ylabel=r'# markers', hist=1, ylim=[0, max_mark], ax=ax)

    def plot_histo_initial_Emax(self, dE=0.5, max_mark=500):
        """Plot initial R of lost particles
        """
        try:
            self.R_w.mean()
        except:
            self._readwall() 

        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        num_wall=np.size(ind)
        x_final = self.data_e['energy'][ind]*1e-3
        x_initial = self.data_i['energy'][ind]*1e-3
        deltaE = x_initial-x_final
        ind_minE = np.where(deltaE<dE)
        num_OL= float(np.size(ind_minE))
        self.plot_histo_E(max_mark=max_mark)
        ax=plt.gca()
        ax.axvline(dE, color='r', lw=2.3, linestyle='--')
        ax.text(-3, 0.86*max_mark, 'FOL \n{:2.1f}%'.format(num_OL*100./num_wall), color='r', fontsize=20)
        ax.text(4, 0.86*max_mark, 'DOL \n', color='b', fontsize=20)

        x = self.data_i['R'][ind[ind_minE]]
        #x = self.data_i['z'][ind[ind_minE]]
        self.plot_histo_R(max_mark=max_mark)
        ax=plt.gca()
        _plot_1d(x, hist=1, ax=ax, color='r')
        ax.axvline(min(self.R_w),color='k', lw=2.3)
        ax.axvline(max(self.R_w),color='k', lw=2.3)

        legend_elements = [patches.Patch(facecolor='b', edgecolor='b',
                                         label=r'$\Delta E>{:.1f} keV$ | {:.1f}%'.format(dE, (1.-num_OL/num_wall)*100.)),

                           patches.Patch(facecolor='r', edgecolor='r',
                                         label=r'$\Delta E<{:.1f} keV$ | {:.1f}%'.format(dE, num_OL/num_wall*100.))]

        ax.legend(handles=legend_elements, loc='upper right')

    def plot_Etime(self, ax=0):
        """
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        x = self.data_e['energy'][ind]
        y = self.data_e['time'][ind]
        
        _plot_2d(x,y, scatter=1, xlabel=r'E [kev]', ylabel=r't [s]', ax=ax)


    def plot_maxrho_histo(self):
        """
        Histogram with final theta position, pitch, energy and 2D plot of the particle velocity
        """
        try:
            self.R_w.mean()
        except:
            self._readwall()
        #ind = np.where(self.data_e['endcond']== 10)[0] #rho=1
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        #ind = np.arange(len(self.data_e['endcond'])) #all particles
        
        pitchi = self.data_i['pitch'][ind]
        energyi = self.data_i['energy'][ind]
        #pitch = self.data_e['pitch'][ind]
        vr = self.data_e['vR'][ind]
        vz = self.data_e['vz'][ind]
        #vphi = self.data_e['vphi'][ind]
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        R0=self.infile['misc/geomCentr_rz'][0]
        theta = np.arctan2(z,r-R0)
        phi = self.data_e['phi'][ind]
        x = r*np.cos(phi); y=r*np.sin(phi)        
      
        
        
        #plt.close('all')
        plt.figure(); plt.hist(pitchi, bins=20); plt.xlabel('Pitch'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(energyi, bins=30); plt.xlabel('Energy'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(vr*1e-3, bins=20); plt.xlabel(r'$v_r$ [km/s]'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(vz*1e-3, bins=20); plt.xlabel(r'$v_z$ [km/s]'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(phi, bins=20); plt.xlabel('Phi (toroidal angle)'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(theta, bins=20); plt.xlabel('theta (poloidal angle)'); plt.ylabel('Number of particles')

        plt.figure(); plt.scatter(x,y);  plt.grid('on'); plt.xlabel(r'x'); plt.ylabel('y')
        plt.tight_layout()

    def save_phitheta_losses(self):
        """
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        R0=self.infile['misc/geomCentr_rz'][0]
        theta = np.arctan2(z,r-R0)
        phi = self.data_e['phi'][ind]        
        hist = np.histogram2d(theta, phi, bins=30)
        f_name_hist='hist_phithetaloss_'+self.id+'.dat'
        with open(f_name_hist, 'w+') as f_handle:
            f_handle.write('phi '+ str(len(hist[1]))+'\n')
            f_handle.write('theta '+str(len(hist[2]))+'\n')
            np.savetxt(f_handle,hist[1]) #phi
            np.savetxt(f_handle,hist[2]) #theta
            np.savetxt(f_handle,hist[0])
        
    def _power_coupled(self):
        """
        Calculates the power coupled to the plasma:
        Power_ini - power_end + power_residual
        """
        p_ini = np.dot(self.data_i['weight'], self.data_i['energy'])*1.602e-19
        p_end = np.dot(self.data_e['weight'], self.data_e['energy'])*1.602e-19
        # you need to add the power which you assume absorbed, i.e. the thermalized, minimum energy and cputime
        endcond = self.data_e['endcond'] 
        flags = [1,2,4]
        ix = np.isin(endcond, flags)
        ind = np.where(ix)[0]
        #ind = np.where(np.logical_or(endcond == 1, endcond == 2, endcond == 4))  #Tmax, Emin
        p_res = np.dot(self.data_e['weight'][ind],self.data_e['energy'][ind])*1.602e-19
        self.pcoup = p_ini-p_end+p_res  
        self.pini = p_ini
        ind = np.where(endcond==3)[0]        
        self.p_lost = np.dot(self.data_e['weight'][ind], self.data_e['energy'][ind])*1.602e-19

    def print_power_coupled(self):
        """
        """
        try:
            self.pcoup
        except:
            self._power_coupled()
            
        print("Injected power ", self.pini*1e-6, " MW")
        print("Coupled  power ", self.pcoup*1e-6, " MW, ", self.pcoup/self.pini*100.," % of Pinj")
         

    def plot_rhopitch(self):
        """ plot rhopitch of losses
    
        Method to plot rho, pitch histo for losses
        
        Parameters:
            None
    
        Returns:
            None
        """
        ind = np.where(self.data_e['endcond']==3)[0]
        rho = self.data_i['rho'][ind]
        pitch = self.data_i['pitch'][ind]
        xlabel = r'$\rho$'; ylabel=r'$\xi=v_\parallel/v$'
        _plot_2d(rho, pitch, hist=1, ylim=[-1, 1], xlabel=xlabel, ylabel=ylabel)
        
        
class SA_iniend(h5_particles):
    def __init__(self, infile_n, fname_surf=''):
        self.device = 'JT60SA'
        self.id = infile_n[-9:-3]
        self.xlim=[1.5, 4.5]
        h5_particles.__init__(self, infile_n, fname_surf)
        
    def plot_beams_XY(self):
        """
        Method to plot XY of ionisation FOR EACH BEAM
        """
        try:
            self.beamorigindict
        except:
            self._calc_originbeam()
            

        #x_range=[-5,5]
        axxy = plt.subplot2grid((1,2),(0,0))
        axrz = plt.subplot2grid((1,2),(0,1))
        ind = np.zeros((2,self.npart), dtype=bool)
        
        for i,el in enumerate(self.beamorigindict):
            if el != 'NNBI_U' and el!='NNBI_L':
                ind[0,:] = self.data_i['origin']==self.beamorigindict[el][0]
                ind[1,:] = self.data_i['origin']==self.beamorigindict[el][1]
            elif el=='NNBI_L':
                continue
            else:
                ind[0,:] = self.data_i['origin']==self.beamorigindict[el]
                ind[1,:] = self.data_i['origin']==self.beamorigindict['NNBI_L']


            for jj in (0,1):
                R   = self.data_i['Rprt'][ind[jj,:]]
                ang = self.data_i['phiprt'][ind[jj,:]]
                z   = self.data_i['zprt'][ind[jj,:]]
                x = np.zeros(len(R))
                y = np.zeros(len(R))
                for j, el in enumerate(R):
                    x[j] = el*math.cos(ang[j])
                    y[j] = el*math.sin(ang[j])
                axxy.scatter(x,y,c=col[i])            
                axrz.scatter(R,z,c=col[i])
                
        theta=np.arange(0,6.29,0.02*6.28)
        if len(self.R_w)==0:
            self._readwall()
        axxy.plot(np.min(self.R_w)*np.cos(theta) , np.min(self.R_w)*np.sin(theta), 'm')
        axxy.plot(np.max(self.R_w)*np.cos(theta) , np.max(self.R_w)*np.sin(theta), 'm')
        axrz.plot(self.R_w , self.z_w, 'm')
        axxy.axis([-5,5,-5,5])
        axrz.axis([0,5,-4,4])
        axxy.axis('equal')
        axrz.axis('equal')

        plt.show()


    def calc_shinethr(self):
        """ Calculate shinethrough

        Method to calculate the shine-through with the weights

        Parameters:
            None
        Returns:
            None        
        """
        
        if self._endgroup != 'shinethr':
            print("WARNING! Check input file, which is ", self._fname)
            
        self._calc_weight()
        id2beamnum = \
                          {\
                           '45':1 ,  '46':1,    '47':2,    '48':2,   \
                           '133':3,  '134':3,   '135':4,   '136':4,  \
                           '221':5,  '222':5,   '223':6,   '224':6,  \
                           '309':7,  '310':7,   '311':8,   '312':8,  \
                           '3637':9, '3638':9,  '3639':10, '3640':10,\
                           '5253':13,'5254':13, '5255':14, '5256':14,\
                           '3031':99,'3032':101 \
                       }
        shinethr = dict.fromkeys(id2beamnum.keys())
        shinethr_abs = dict.fromkeys(id2beamnum.keys())

        for i in id2beamnum:
            ind = self.data_e['origin'][:]==int(i)
            if len(self.data_e['origin'][ind])==0:
                w=0
                e=0
                power = 1
            else:  
                power=1.e6
                if int(i) in [3031, 3032]:
                    power = 5.e6
            try:
                e = self.data_e['energy'][ind][0]
            except:
                e=0
            e *= 1.612e-19
            w = np.sum(self.data_e['weight'][ind])
            #wtot = self.w[ind_i]
            shinethr[i]=float(e)*w/power
            shinethr_abs[i] = float(e)*w

        beamnum2id = \
                     {'1':[45, 46],      '2':[47, 48],\
                      '3':[133, 134],    '4':[135, 136],\
                      '5':[221, 222],    '6':[223, 224],\
                      '7':[309, 310],    '8':[311, 312],\
                      '9':[3637, 3638],  '10':[3639, 3640],\
                      '13':[5253, 5254], '14':[5255, 5256],\
                      '99':[3031],       '101':[3032]}  
        self.shinethr     = dict.fromkeys(beamnum2id.keys(), 0)
        self.shinethr_abs = dict.fromkeys(beamnum2id.keys(), 0)
        for i in beamnum2id:
            for j in range(len(beamnum2id[i])):
                self.shinethr[i]+=shinethr[str(beamnum2id[i][j])]
                self.shinethr_abs[i]+=shinethr_abs[str(beamnum2id[i][j])]
                
class TCV_iniend(h5_particles):
    def __init__(self, infile_n):
        self.device = 'TCV'
        self.id = infile_n[-11:-5]
        self.xlim = [0.6, 1.1]
        h5_particles.__init__(self, infile_n)
		
    def calc_shinethr(self):
         """
         Method to calculate the shine-through with the weights
         """

         if self._endgroup != 'shinethr':
             print("WARNING! Check input file, which is ", self._fname)

         self._calc_weight()
         self.shinethr=np.zeros((1), dtype=float)
         self.shinethr_abs=np.zeros((1), dtype=float)
         power = 0.62e6

         e = self.data_e['energy']
         e = e*1.612e-19
         w = self.data_e['weight']
         self.shinethr_abs= np.sum(e*w)
         self.shinethr=np.sum(e*w)/power

         print("TCV Shine-through:", "{0:5f} %".format(self.shinethr*100),\
                                   "||  {0:3f} W".format(self.shinethr_abs))

def _readwall(device):
    """
    Hidden method to read the wall
    """
    if device == 'JT60SA':
        in_w_fname = '/home/vallar/ASCOT/runs/JT60SA/002/input.wall_2d'
    elif device == 'TCV':
        cluster = platform.uname()[1]
        if cluster[-7:] == 'epfl.ch':
            in_w_fname = '/home/vallar/TCV/input.wall_2d_FW'
        else:
            in_w_fname = '/home/vallar/ASCOT/runs/TCV/57850/input.wall_2d'

    wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=1)
            
    R_w = wall[0,:]
    z_w = wall[1,:]
    R_w = np.array(R_w)
    z_w = np.array(z_w)
    return R_w, z_w
		
def _RZsurf(fname_surf=''):
    """
    Reads the position of RZ surfaces from ascot file
    """               
    if fname_surf=='':
        f = h5py.File('ascot.h5')
        fname_surf = ' ascot.h5 '
    else:
        f = h5py.File(fname_surf)

    print("READING SURF FROM "+fname_surf)#ascot_"+str(id)+".h5")
    RZsurf = -1.*f['bfield/2d/psi'].value
    Rsurf = f['bfield/r']
    zsurf = f['bfield/z']
    try:
        edge = f['boozer/psiSepa'][:]; axis=f['boozer/psiAxis'][:]
    except:
        edge=1; axis=0 
    print('Psi edge : ', edge, ' Psi axis : ', axis)
    RZsurf = (RZsurf-axis)/(edge-axis)
    RZsurf = np.sqrt(RZsurf)   
    return Rsurf, zsurf, RZsurf
	

