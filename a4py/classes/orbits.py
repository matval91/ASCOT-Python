from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, gridspec
import collections
import math, platform
import h5py, random
from mpl_toolkits.mplot3d import Axes3D
from utils.plot_utils import _plot_2d, _plot_RZsurf, limit_labels



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

class dat_particles(particles):
    """
    superClass (inherited from particles) handling dat files (e.g. input.particles, test_ascot.orbits.dat)
    """
    def __init__(self, infile_n, fname_surf=''):
        """
        READ ASCII FILE
        """
        self._fname = infile_n 
        self.id = infile_n[-10:-4]

        if fname_surf=='':
            fname_surf = 'ascot_'+self.id+'.h5'
        self.infile = h5py.File(fname_surf)
        particles.__init__(self, fname_surf)
        in_f  = open(self._fname)
        lines = in_f.readlines()[3:]
        in_f.close()
        #Read lines of comment
        n_comm = int(lines[0].split()[0])
        
        # read number of particles
        self.npart = int(lines[n_comm+2].split()[0])

        # Number of fields
        self.nfields = int(lines[n_comm+4].split()[0])
        line_fields  = n_comm+5
        # Read and store the fields and their unit
        for i in range(self.nfields):
            self.field.append(lines[line_fields+i].split()[0])
            self.unit.append(lines[line_fields+i].split()[-1][1:-1])
        self.field = np.array(self.field)
        self.unit = np.array(self.unit)
        if self.npart==-1:
            self.npart=input("Number of particles? ")
            self.npart=int(self.npart)
        tmpdict = dict.fromkeys(self.field,[])
        self.partdict = np.array([dict.copy(tmpdict) for x in range(self.npart)])
        
        # Store actual data 
        part = lines[line_fields+self.nfields+1:-1]
        ind_idfield = np.argwhere(self.field == 'id')[0]
        for rowt in part:
            row = rowt.split()
            row = np.array(row, dtype=float)
            part_id = int(row[ind_idfield][0]-1)
            if part_id%5000==0:
                print("Doing particle ", part_id)
            for i, el in enumerate(row):
                if self.field[i]=='phi':
                    self.partdict[part_id][self.field[i]] = \
                    np.append(self.partdict[part_id][self.field[i]], el*math.pi/180.)                   
                else:
                    self.partdict[part_id][self.field[i]] = \
                    np.append(self.partdict[part_id][self.field[i]], el)                

        self.data_i = dict.copy(tmpdict)
        self.data_e = dict.copy(tmpdict)

        for i in range(self.npart):
            for key in self.data_i:
                self.data_i[key] = np.append(self.data_i[key], \
                            self.partdict[i][key][0])
                self.data_e[key] = np.append(self.data_e[key], \
                                             self.partdict[i][key][-1])
#        for key in self.field:
#            if key=='id':
#                self.data_i[key] = np.array(self.data_i[key], dtype=int)
#            else:
#                self.data_i[key] = np.array(self.data_i[key], dtype=float)

        #CONVERT PHI FROM DEG TO RAD
        #self.data_i['phi'] *= math.pi/180.
        #self.origins = np.array(list(set(self.data_i['origin'])))
        print("FIELDS IN FILE ",self._fname," :")
        print(self.field)
        self.banana_orbit()
         
    def plot_histo_wall(self, ax=0, ini=0):
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
        ylim = [-2.4, -0.8]
        #ylim = [-3.14, 3.14]
        _plot_2d(phi, theta, xlabel=r'$\phi$ [rad]',ylabel=r'$\theta$ [rad]', Id = self.id, hist=1, xlim=[-3.14, 3.14],ylim=ylim)
        ax=plt.gca()

        if ini!=0:
            r = self.data_i['R'][ind]
            z = self.data_i['z'][ind]
            theta = np.arctan2(z-z0,r-R0)
            phi = self.data_i['phi'][ind]
            _plot_2d(phi, theta, xlabel='', ylabel='', ax=ax, scatter=1)
    
    def plot_rhopitch(self, ax=0):
        """
        Histogram of deposition to the wall
        """
        ind = np.where(self.trappind==1)[0]
        rho = self.data_i['rho'][ind]
        pitch = self.data_i['pitch'][ind]
        
        _plot_2d(rho, pitch, scatter=1, xlabel=r'$\rho$', ylabel=r'$\xi$')

    def _compute_npart(self, line):
        """
        computes the number of particles if it is not possible to gather it from
        the dat file. Takes the last line and checks the particle's id
        """
        line = line.split()
        part_id = int(line[np.argwhere(self.field == 'id')[0]])
        self.npart = part_id
                
    def banana_orbit(self):
        """
        Computes the fraction of the particles doing a trapped orbit
        To do it, set the simulation time to 1e-3, 1e-4 and the writeInterval small.
        from the orbits, you check where the pitch changes sign: this means the orbit is trapped
        trappind=1 means orbit is trapped
        """         
        self.ntrapp = 0
        self.trappind = np.zeros(self.npart)
        for i in range(self.npart):
            pitch=self.partdict[i]['pitch'][:]
            pitch_s = np.sign(pitch)
            signchange = ((np.roll(pitch_s, 1) - pitch_s) != 0).astype(int)
            signchange = np.array(signchange)
            if len(np.where(signchange==1)[0])!=0:
                self.ntrapp+=1
                self.trappind[i]=1
                if i<0:
                    print(len(np.where(signchange==1)[0]))
                    f=plt.figure(); f.suptitle(i)
                    ax=f.add_subplot(211); ax.plot(self.partdict[i]['R'], self.partdict[i]['z'], 'k-')
                    ax2=f.add_subplot(212); ax2.plot(pitch);ax2.plot(signchange)
                    plt.ginput()
                    plt.close('all')
                
        self._banana_dim()

    def plot_random_orbit(self, npart=20, ind=[0]):
        """ Plots orbit of some random particles
        
        Initialises some random particles and plots their orbits, just to check what's happening
        
        Parameters:
            None
        Returns:
            None
        """
        f, axtrajxy, axtrajRZ = self._initialize_figure_ptcls_2d()
        f3d, ax3d = self._initialize_figure_ptcls_3d()
        if npart==20:
            for ii in range(npart):
                ind2plot = random.randint(0, self.npart)
                p2plot = self.partdict[ind2plot]
                #actual plot
                self._plot_trajectory(p2plot=p2plot, f=f, axRZ=axtrajRZ, axxy=axtrajxy)
                self._plot_3d_trajectory(p2plot=p2plot, ax=ax3d)
            axtrajxy.set_xlim([-max(self.R_w), max(self.R_w)])
            axtrajRZ.set_xlim([min(self.R_w)*0.9, max(self.R_w)*1.1])
        else:
            col=['k', 'r','g','b','m','c']
            for ii in range(len(ind)):
                ind2plot = ind[ii]
                p2plot = self.partdict[ind2plot]
                #actual plot
                self._plot_trajectory(p2plot=p2plot, f=f, axRZ=axtrajRZ, axxy=axtrajxy, col=col[ii])
                self._plot_3d_trajectory(p2plot=p2plot, ax=ax3d,col=col[ii])
            axtrajxy.set_xlim([-max(self.R_w), max(self.R_w)])
            axtrajRZ.set_xlim([min(self.R_w)*0.9, max(self.R_w)*1.1])  
          
    def plot_orbit_wall(self, theta_f=0., phi_f=0., npart=20):
        """
        Chosen a position on the wall, gets the first ten particles which end there and plots their orbit
        """
        #Filtering on the range desired in phi, theta
        if theta_f==0. and phi_f==0.:
            self.plot_histo_wall()
            ax_histo=plt.gca()
            p1, p2 = plt.ginput(n=2)
            theta_f = [p1[1], p2[1]]
            phi_f = [p1[0], p2[0]] 

        phi_f   = np.linspace(phi_f[0], phi_f[1], npart)
        theta_f = np.linspace(theta_f[0], theta_f[1], npart)

        #==========================================================
        # Finding index where particle meet the conditions
        #==========================================================        
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]        
        R0 = self.infile['misc/geomCentr_rz'][0]
        z0 = self.infile['misc/geomCentr_rz'][1]
        theta = np.arctan2(z-z0,r-R0)
        phi = self.data_e['phi'][ind]
        ax_histo.scatter(phi_f, theta_f, 100, marker='*', c='r')
        Rw=self.R_w
        zw=self.z_w

        f, axtrajxy, axtrajRZ = self._initialize_figure_ptcls_2d()
        f3d, ax3d = self._initialize_figure_ptcls_3d()
        f_rhopitch = plt.figure(); axrhopitch=f_rhopitch.add_subplot(111)

        for ii in range(len(phi_f)):
            # finding particles that are closes to the points chosen
            dist = np.sqrt((theta-theta_f[ii])**2+(phi-phi_f[ii])**2)
            ind_inist = np.argmin(dist)
            ind_inistate = ind[ind_inist]
            
            #index of particles to plot
            p2plot = self.partdict[ind_inistate]
            
            #plotting initial phi,theta in histogram
            R_ini,z_ini, phi_ini = p2plot['R'][0], p2plot['z'][0], p2plot['phi'][0]
            theta_ini = np.arctan2(z_ini-z0,R_ini-R0)
            ax_histo.scatter(phi_ini, theta_ini, s=20, marker='+', color='r')
            axrhopitch.scatter(p2plot['rho'][0], p2plot['pitch'][0], s=20)
            #actual plot
            self._plot_trajectory(p2plot=p2plot, f=f,axRZ=axtrajRZ, axxy=axtrajxy)
            self._plot_3d_trajectory(p2plot=p2plot, ax=ax3d)
        
        axtrajxy.set_xlim([-max(self.R_w), max(self.R_w)])
        axtrajRZ.set_xlim([min(self.R_w), max(self.R_w)])

    def _banana_dim(self):
        """
        computes the fraction of trapped particles with the following formula:
         ft = sqrt(epsilon)*(1.46-0.46*epsilon)
         
        and the width of the banana with the following formula
         w_b = sqrt(epsilon)*rho_L(larmor radius evaluated using only poloidal field)
        w_b =2 v_para m/ (q B_p[(a+r)/2])
        """
        try:
            self.partdict[0]['mass'].mean()
        except:
            print("Impossible to calculate the banana orbits dimension")
            return
        try:
            self.partdict[0]['charge'].mean()
        except:
            for i in self.partdict:
                i['charge']=1
            
        m, E = self.data_i['mass'], self.data_i['energy']
        q = np.full(len(m),1.)*1.602e-19; E = E*1.602e-19
        m = m*1.66e-27
        v = np.sqrt(2.*E/m)*self.data_i['pitch']
        r = self.data_i['R']-self.R0
        Bpol = np.sqrt(self.data_i['BR']**2+self.data_i['Bz']**2)
        #print(Bpol.shape)
        w_banana = np.abs(v)*m/(q*Bpol)
        self.w_banana_all = w_banana[np.where(self.trappind==1)[0]]
        self.w_banana = w_banana[np.where(np.logical_and(self.trappind==1, self.data_i['rho']>0.4))[0]]


    def plot_trapped_contour(self):
        """
        plots trapped particles in different spaces: (R,z), (vpara, vperp),
        (pitch, energy)
        RED: trapped
        BLUE: not trapped
        """
        try:
            self.ntrapp
        except:
            self.banana_orbit()
            
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)    
        #plt.rc('font', size=20)

        ylim=[-1., -0.5]
        
        ind_t = np.where(self.trappind==1)[0]
        ind_p = np.where(self.trappind==0)[0]     
        CL = ['','']
        r = self.data_i['R']; z = self.data_i['z']
        p = self.data_i['pitch']
        rho = self.data_i['rho']
        qpart = np.zeros(self.npart)
            
        nbins=10; nsurf=20
        hist_p, xedges_p, yedges_p = np.histogram2d(r[ind_p], z[ind_p], \
                                    bins=nbins)
        hist_t, xedges_t, yedges_t = np.histogram2d(r[ind_t], z[ind_t], \
                                    bins=nbins)
        vmaxt = max(np.max(hist_t), np.max(hist_p))
        ind_cb=1
        if np.max(hist_p)> np.max(hist_t):
            ind_cb=0
            
        f=plt.figure(figsize=(11,10))
        #f.suptitle(self.id)
        f.text(0.01, 0.95, self.id)
        
        axrz = f.add_subplot(221)
        x = np.linspace(np.min(r[ind_p]), np.max(r[ind_p]), num=nbins)
        y = np.linspace(np.min(z[ind_p]), np.max(z[ind_p]), num=nbins)
        CL[0]=axrz.contourf(x,y,hist_p.T, nsurf, cmap=my_cmap, vmin=0, vmax=vmaxt)
#        CL[0]=axrz.pcolor(x,y,hist_p.T, cmap=my_cmap, vmin=0, vmax=vmaxt)
        
        _plot_RZsurf(self.Rsurf, self.zsurf, self.RZsurf,axrz)
        axrz.set_title('Passing N(R,Z)'); 
        axrz.set_xlabel('R [m]'); axrz.set_ylabel('Z [m]')        
        axrz.axis('equal');         
#        axrz.set_xlim([np.min(x), np.max(x)]);
        axrz.set_ylim([-1., 1.])
        axrz2 = f.add_subplot(222)
        x = np.linspace(np.min(r[ind_t]), np.max(r[ind_t]), num=nbins)
        y = np.linspace(np.min(z[ind_t]), np.max(z[ind_t]), num=nbins)
        CL[1]=axrz2.contourf(x,y,hist_t.T, nsurf, cmap=my_cmap, vmin=0, vmax=vmaxt)
#        CL[1]=axrz2.pcolor(x,y,hist_t.T, cmap=my_cmap, vmin=0, vmax=vmaxt)

        _plot_RZsurf(self.Rsurf, self.zsurf, self.RZsurf,axrz2)
        axrz2.set_title('Trapped N(R,Z)'); 
        axrz2.set_xlabel('R [m]'); axrz2.set_ylabel('Z [m]')
        cbar_ax = f.add_axes([0.85, 0.6, 0.03, 0.3])
        f.colorbar(CL[ind_cb], cax=cbar_ax)
        axrz2.axis('equal')
#        axrz2.set_xlim([np.min(x), np.max(x)]); 
#        axrz2.set_ylim([np.min(y), np.max(y)])
        axrz2.set_ylim([-1., 1.])


        hist_t, yedges_t, xedges_t = np.histogram2d(rho[ind_t], p[ind_t],bins=nbins)
        hist_p, yedges_p, xedges_p = np.histogram2d(rho[ind_p], p[ind_p],bins=nbins)
        vmaxt = max(np.max(hist_t), np.max(hist_p))   

        f2=plt.figure()
        f2.suptitle(self.id)
        ax2=f2.add_subplot(111)        
        axrp = f.add_subplot(223)
        x = np.linspace(np.min(rho[ind_p]), np.max(rho[ind_p]), num=nbins)
        y = np.linspace(np.min(p[ind_p]), np.max(p[ind_p]), num=nbins)
        CL[0]=axrp.contourf(x,y,hist_p.T, nsurf, cmap=my_cmap, vmin=0, vmax=vmaxt)
#        CL[0]=axrp.pcolor(x,y,hist_p.T, cmap=my_cmap, vmin=0, vmax=vmaxt)

        ax2.contour(x,y,hist_p.T,nsurf,colors='k', label='Passing')
        
        axrp.set_title(r'Passing N($\rho$, $\xi$)'); axrp.set_xlabel(r'$\rho$'); axrp.set_ylabel(r'$\xi$')  
#        axrp.set_xlim([np.min(x), np.max(x)]); axrp.set_ylim([np.min(y), np.max(y)])
        axrp.set_xlim([0., 1.]); axrp.set_ylim(ylim)
        
        axrp2 = f.add_subplot(224)
        x = np.linspace(np.min(rho[ind_t]), np.max(rho[ind_t]), num=nbins)
        y = np.linspace(np.min(p[ind_t]), np.max(p[ind_t]), num=nbins)
        CL[1] = axrp2.contourf(x,y,hist_t.T, 10, cmap=my_cmap, vmin=0, vmax=vmaxt)
#        CL[1] = axrp2.pcolor(x,y,hist_t.T, cmap=my_cmap, vmin=0, vmax=vmaxt)
        ax2.contour(x,y,hist_t.T,nsurf,colors='r', label='Trapped')
        axrp2.set_title(r'Trapped N($\rho$, $\xi$)'); axrp2.set_xlabel(r'$\rho$'); axrp2.set_ylabel(r'$\xi$')  
#        axrp2.set_xlim([np.min(x), np.max(x)]); axrp2.set_ylim([np.min(y), np.max(y)])  
        axrp2.set_xlim([0., 1.]); axrp2.set_ylim(ylim)      
        
        ax2.set_xlabel(r'$\rho$'); ax2.set_ylabel(r'$\xi$')
        ax2.legend(loc='upper right')
        for aa in [axrz, axrz2, axrp, axrp2]:        
            aa.xaxis.set_major_locator(plt.MaxNLocator(4))
            aa.yaxis.set_major_locator(plt.MaxNLocator(4))            
        f.tight_layout()
        f.subplots_adjust(right=0.8)
        cbar_ax = f.add_axes([0.85, 0.1, 0.03, 0.35])
        f.colorbar(CL[ind_cb], cax=cbar_ax)
        axrp.grid('on')
        axrp2.grid('on')

    def _initialize_figure_ptcls_2d(self):
        """ initialize figure for more particles
        
        Private method to initialize the figure in 2d and 3D for orbit plot

        Parameters:
            None
        Returns:
            | f
            | axtrajxy
            | axtrajRZ
            | f3d
            | ax3d
        
        """
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        plt.rc('figure', facecolor='white')
        #defining figure
        f=plt.figure(figsize=(12,9))
        gs = gridspec.GridSpec(1,2, width_ratios=[1,3])
        axtrajRZ = f.add_subplot(gs[0])
        axtrajxy = f.add_subplot(gs[1])
        #plots RZ wall and surfaces
        _plot_RZsurf(self.Rsurf, self.zsurf, self.RZsurf, axtrajRZ)
        axtrajRZ.plot(self.R_w, self.z_w, 'k', lw=2.5)
        # plots magnetic axis
        circle1 = plt.Circle((0, 0), self.R0, color='r', fill=False, linestyle='--', linewidth=3.)      
        axtrajxy.add_artist(circle1)
        #plots wall in XY
        circle1 = plt.Circle((0, 0), min(self.R_w), color='k', fill=False, linestyle='-', linewidth=3.)      
        axtrajxy.add_artist(circle1)
        circle1 = plt.Circle((0, 0), max(self.R_w), color='k', fill=False, linestyle='-', linewidth=3.)      
        axtrajxy.add_artist(circle1)
        axtrajxy.set_xlabel(r'x (m)'); axtrajxy.set_ylabel(r'y (m)')     
        axtrajRZ.axis('equal'); axtrajxy.axis('equal')
        axtrajxy.set_xlim([-max(self.R_w), max(self.R_w)])
        axtrajRZ.set_xlim([min(self.R_w), max(self.R_w)])
        limit_labels(axtrajRZ, r'R (m)', r'Z (m)', M=4)
        axtrajxy.grid('on')
        return f, axtrajxy, axtrajRZ
        
    def _initialize_figure_ptcls_3d(self):

        f3d=plt.figure(figsize=(15,9))
        ax3d = f3d.add_subplot(111, projection='3d')
        #plot 3d tokamak 
        phi_tt = np.arange(0.*np.pi,2.02*np.pi,0.02*np.pi)
        the_tt = phi_tt
        x_tok = np.zeros((len(phi_tt),len(self.R_w)),dtype=float)
        y_tok = np.zeros((len(phi_tt),len(self.R_w)),dtype=float)
        for i,R in enumerate(self.R_w):
            x_tok[:,i] = R*np.cos(phi_tt)
            y_tok[:,i] = R*np.sin(phi_tt)
        z_tok = np.array(self.z_w)
        z_tok=np.tile(z_tok, [np.size(phi_tt),1])
        ax3d.plot_surface(x_tok,y_tok,z_tok,color='k',alpha=0.15)

        return f3d, ax3d

        
    def _plot_trajectory(self, p2plot=0, f=0, axRZ=0, axxy=0, col='k'):
        """
        Plots trajectory of one single particle
        """
        try:
            self.trappind.mean()
        except:
            self.banana_orbit()

        #if a particles is desired, it choses that particle
        if p2plot == 0.:
            ind2plot = 0
        else:
            ind2plot = p2plot['id'][0]-1
        print('plotting particle '+str(ind2plot+1))

        if f==0:
            f, axtrajxy, axtrajRZ = self._initialize_figure_ptcls_2d()
        else:
            axtrajRZ = axRZ
            axtrajxy = axxy
        ind2plot=int(ind2plot)
        self._plot_traj_RZ(axtrajRZ, ind2plot, col=col)
        self._plot_traj_xy(axtrajxy, ind2plot, col=col)
        f.tight_layout()
        plt.show()

    def _plot_traj_RZ(self, axtrajRZ, ind2plot, col='k'):
        #plots RZ ini and enstate
        axtrajRZ.scatter(self.data_i['R'][ind2plot], self.data_i['z'][ind2plot], 100, 'k', marker='*')
        axtrajRZ.scatter(self.data_e['R'][ind2plot], self.data_e['z'][ind2plot], 100,'r', marker='*')

        # plot RZ trajectory
        axtrajRZ.plot(self.partdict[ind2plot]['R'], self.partdict[ind2plot]['z'], c=col, linestyle='-', lw=2.3)    


    def _plot_traj_xy(self, axtrajxy, ind2plot, col='k'):
        #calculates xy trajectory
        xi,yi = self.data_i['R'][ind2plot]*math.cos(self.data_i['phi'][ind2plot]), \
                self.data_i['R'][ind2plot]*math.sin(self.data_i['phi'][ind2plot])
        x,y = self.partdict[ind2plot]['R']*np.cos(self.partdict[ind2plot]['phi']), \
              self.partdict[ind2plot]['R']*np.sin(self.partdict[ind2plot]['phi'])  
        xe,ye = self.data_e['R'][ind2plot]*np.cos(self.data_e['phi'][ind2plot]), \
                self.data_e['R'][ind2plot]*np.sin(self.data_e['phi'][ind2plot])                  #plots xy ini/endstate and trajectory
        axtrajxy.scatter(xi,yi, 100, 'k',marker='*')
        axtrajxy.plot(x,y, c=col, linestyle='-', lw=2.3)
        axtrajxy.scatter(xe, ye, 100,'r', marker='*')        

    def _plot_3d_trajectory(self, p2plot=0, ax=0, col='k'):
        """
        Plots trajectory of one single particle in 3D
        """
        try:
            self.trappind.mean()
        except:
            self.banana_orbit()

        #if a particles is desired, it choses that particle
        if p2plot == 0.:
            ind2plot = 0
        else:
            ind2plot = p2plot['id'][0]-1
        ind2plot=int(ind2plot)
        print('plotting particle '+str(ind2plot+1))
        
        if ax==0:
            f = plt.figure(figsize=(10,5))
            f.suptitle('Orbit particle '+str(ind2plot+1))
            ax3d = fig3d.add_subplot(111, projection='3d')
            ##3D plot
            #TOKAMAK
            #shape of tokamak
            phi_tt = np.arange(0.*np.pi,2.02*np.pi,0.02*np.pi)
            the = phi
            x_tok = np.zeros((len(phi),len(self.R_w)),dtype=float)
            y_tok = np.zeros((len(phi),len(self.R_w)),dtype=float)
            for i,R in enumerate(self.R_w):
                x_tok[:,i] = R*np.cos(phi)
                y_tok[:,i] = R*np.sin(phi)
            z_tok = self.z_w
            ax3d.plot_surface(x_tok,y_tok,z_tok,color='k',alpha=0.15)

        else:
            ax3d=ax

        #plots RZ ini 
        R = self.data_i['R'][ind2plot]
        z = self.data_i['z'][ind2plot]
        phi = self.data_i['phi'][ind2plot]
        x, y = R*np.cos(phi), R*np.sin(phi)
        ax3d.scatter(x,y, zs=z, s=20, c='k', marker='*')

        #plots RZ end
        R = self.data_e['R'][ind2plot]
        z = self.data_e['z'][ind2plot]
        phi = self.data_e['phi'][ind2plot]
        x, y = R*np.cos(phi), R*np.sin(phi)
        ax3d.scatter(x,y, zs=z, s=20, c='k', marker='*')

        # plot RZ trajectory
        R,z = self.partdict[ind2plot]['R'], self.partdict[ind2plot]['z']
        phi = self.partdict[ind2plot]['phi']
        x, y = R*np.cos(phi), R*np.sin(phi)
        ax3d.plot(x,y, zs=z, c=col, linestyle='-', lw=2.5)

        ax3d.set_xticklabels([''])
        ax3d.set_yticklabels([''])
        ax3d.set_zticklabels([''])

        ax3d.set_xlabel('X'); ax3d.set_ylabel('Y');
        ax3d.set_zlabel('Z');
        plt.axis('off')
        plt.show()


    def detrapping(self):
        """
        Calculates the detrapping condition for the particles
        As calculated in the Wesson book (3.12.11), the condition for the detrapping
        to happen is
            tau_coll \leq (R_0/r)^{1.5} q R_0/(sqrt(2)v_perp)
        where
            v_perp = sqrt(2)* v_thermal = sqrt(2)*sqrt(kB T/m)
        The tau_coll to use for comparison is the spitzer E
        """
        try:
            self.taucoll_e.mean()
        except:
            self._colltimes()
        kB = 1.38064852e-23
        me = 9.10938356e-31
        mp = 1.6726219e-27
        e = 1.602e-19
        R_torus = 2.96 #in meters
        a = 1.11
        rho = self.infile['/boozer/psi'][:]**0.5
        q = self.infile['boozer/qprof'][:]*(-1)
        self.epsilon = a/R_torus
        print("Epsilon used is for scenario 5")           
        te = self.infile['plasma/1d/te'][:]      
        vperp = np.zeros(self.npart)   
        self.tau_detrapp = np.zeros(self.npart)

        for i in range(self.npart):        
            E = self.partdict[i]['energy'][0]*e
            m = self.partdict[i]['mass'][0]*mp
            v = math.sqrt(2.*E/m)
            angle = np.arccos(self.partdict[i]['pitch'][0])
            vperp[i] = v*math.sin(angle)
            R,z = self.partdict[i]['R'][0], self.partdict[i]['z'][0]            
            r = math.sqrt((R-R_torus)**2+z**2)
            factor = (R_torus/r)**1.5*R_torus/(math.sqrt(2)*vperp[i])
            ind = np.argmin(rho-self.partdict[i]['rho'][0]>0)
            self.tau_detrapp[i] = factor*q[ind]

    def _ecrit(self):
        """
        Calculates critical energy profiles
        Ec = 
        ts = 6.28e14*(A*te^1.5)/(Z^2*ne*lnlambda)
        """
        rho = self.infile['plasma/1d/rho'][:]
        te = self.infile['plasma/1d/te'][:]
        ne = self.infile['plasma/1d/ne'][:]
        Ai = self.infile['plasma/anum'][:]
        Zi = self.infile['plasma/znum'][:]
        nimp = self.infile['plasma/1d/ni'][:]   
        A = self.infile['species/testParticle/anum'][0]
        Z = self.infile['species/testParticle/znum'][0]
        summ = np.sum(np.multiply(nimp, Zi**2/Ai), axis=1)
        Ec = 14.8*te*(A**(1.5)/ne*summ)**(2./3.)
        self.param_ec = interpolate.interp1d(rho,Ec)

        #Spitzer slowing-down time
        ts = 6.28e14*A*te**1.5/(Z**2*ne*17.)
        self.param_ts = interpolate.interp1d(rho,ts)


class SA_orbits(dat_particles):
    def __init__(self, infile_n, fname_surf=''):
        self.device = 'JT60SA'
        self.id = infile_n[-10:-4]
        dat_particles.__init__(self, infile_n, fname_surf)
        self.a = 1.11
        self.R0 = 2.96

    def plot_trapped_energy_PNBs(self):
        """
        plots trapped particles as function of energy
        """
        plt.rc('font', weight='bold')
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
        plt.rc('figure', facecolor='white')
        plt.rc('legend', fontsize=20)

        e = self.data_i['energy']; 
        width=0.2
        num_t_full = len(np.where(np.logical_and(self.trappind==1, e>80000))[0])
        num_t_half = len(np.where(np.logical_and(self.trappind==1, np.logical_and(e>30000, e<80000)))[0])
        num_t_thir = len(np.where(np.logical_and(self.trappind==1, e<30000))[0])
        num_full = len(np.where(e>80000)[0])
        num_half = len(np.where(np.logical_and(e>30000, e<80000))[0])
        num_thir = len(np.where(e<30000)[0])
        num_tot = num_full+num_half+num_thir
        f = plt.figure()
        #f.suptitle(self.id)
        ax = f.add_subplot(111)
        #f.text(0.01, 0.95, str(float(self.ntrapp)/self.npart))
        x = [85000./3.,85000./2.,85000.]
        x = [0.33-width/2., 0.66-width/2., 1.-width/2.]
        y = [float(num_t_thir)/num_thir, float(num_t_half)/num_half, float(num_t_full)/num_full]
        y1 = [float(num_thir)/num_tot, float(num_half)/num_tot, float(num_full)/num_tot]
        y2 = [float(num_t_thir)/num_tot, float(num_t_half)/num_tot, float(num_t_full)/num_tot]

        ax.bar(x, y1, width, color='b', label='Passing')
        ax.bar(x, y2, width, color='r', label='Trapped')
        ax.set_ylim([0, 1.2])
        ax.set_xticks([0.34, 0.66, 1.])
        ax.set_xticklabels([r'$E_{inj}$/3', r'$E_{inj}$/2', r'$E_{inj}$'])
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.])
        
        for i in range(3):
            ax.text(x[i]-0.01, y1[i]+0.05, '{:.1f} %'.format(y[i]*100.), color='r', fontsize=30)

        ax.legend(loc='best')
        ax.set_ylabel(r'Fraction to total number')
        ax.yaxis.grid('on')
        f.tight_layout()

    def plot_trapped_energy_NNBs(self):
        """
        plots trapped particles as function of energy
        """
        plt.rc('font', weight='bold')
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
        plt.rc('figure', facecolor='white')
        plt.rc('legend', fontsize=20)

        e = self.data_i['energy']; 
        width=0.1
        num_t_full = len(np.where(self.trappind==1)[0])
        num_full = len(e)

        f = plt.figure(figsize=(5,8))
        ax = f.add_subplot(111)
        ax.text(0.01, 0.01, self.id)
        x = [500e3]
        x = [0.5]
        y = [float(num_t_full)/num_full]

        ax.bar(x, [1.], width, color='b', label='Passing')
        ax.bar(x, y, width, color='r', label='Trapped')
        ax.set_ylim([0, 1.2])
        ax.set_xlim([0.45, 0.65])
        ax.legend(loc='best')
        ax.yaxis.grid('on')
        ax.set_xticks([0.55])
        ax.set_xticklabels([r'$E_{inj}$'])
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.])
        ax.set_ylabel(r'Fraction to total number')
        ax.text(x[0]-0.05, 1.05, '{:.1f} %'.format(y[0]*100.), color='r')
        f.tight_layout()
        
    def _colltimes_PNBs(self):
        """
        Calculating collision times between fast ions and e/i,
        as in formula 2.15.3 in Wesson book, putting Ec instead of Te/i        
            tau_e = 1.09e16 * T_e[keV]^1.5 / (n_i Z^2 lnlambda) sec
            tau_i = 6.6e17 * (m_i/m_p)^{0.5} T_i[keV]^1.5 / (n_i Z^4 lnlambda) sec
        where m_p should be the proton mass
        
        """
        
        rho = self.infile['plasma/1d/rho'][:]
        E = np.array([85, 85/2., 85/3.]) #in keV    
#        te = np.trapz(te, rho)
        ni = self.infile['plasma/1d/ni'][:,0]   
#        ni = np.trapz(ni, rho)
        Zi = self.infile['plasma/znum'][0]
        lnlambda=17
        Ai = self.infile['plasma/anum'][0]
#        ti = self.infile['plasma/1d/ti'][:]*1e-3
#        ti = np.trapz(ti, rho)
        taucoll_i = np.zeros((np.shape(E)[0], np.shape(ni)[0]))
        taucoll_e = np.zeros((np.shape(E)[0], np.shape(ni)[0]))
        self.taucoll_e = np.zeros((np.shape(E)[0], len(np.where(rho<1.)[0])))
        self.taucoll_i = np.zeros((np.shape(E)[0], len(np.where(rho<1.)[0])))
        self.taucoll_e_mean = np.zeros((np.shape(E)[0]))
        self.taucoll_i_mean = np.zeros((np.shape(E)[0]))
        
        for en_ind, en in enumerate(E):
            taucoll_e[en_ind,:] = 1.09e16*en**1.5/(ni*Zi**2*lnlambda)
            taucoll_e[en_ind, rho>=1.] = np.NaN
            self.taucoll_e[en_ind,:] = taucoll_e[en_ind,~np.isnan(taucoll_e[en_ind,:])]
            taucoll_i[en_ind,:] = 6.6e17*Ai**0.5*en**1.5/(ni*Zi**4*lnlambda)
            taucoll_i[en_ind, rho>=1.] = np.NaN
            self.taucoll_i[en_ind,:] = taucoll_i[en_ind,~np.isnan(taucoll_i[en_ind,:])]
            self.taucoll_e_mean[en_ind] = np.trapz(self.taucoll_e[en_ind,:], rho[rho<1.])
            self.taucoll_i_mean[en_ind] = np.trapz(self.taucoll_i[en_ind,:], rho[rho<1.])        

    def _colltimes_NNBs(self):
        """
        Calculating collision times between fast ions and e/i,
        as in formula 2.15.3 in Wesson book, putting Ec instead of Te/i        
            tau_e = 1.09e16 * T_e[keV]^1.5 / (n_i Z^2 lnlambda) sec
            tau_i = 6.6e17 * (m_i/m_p)^{0.5} T_i[keV]^1.5 / (n_i Z^4 lnlambda) sec
        where m_p should be the proton mass
        
        """
        
        rho = self.infile['plasma/1d/rho'][:]
        E = np.array([500]) #in keV    
#        te = np.trapz(te, rho)
        ni = self.infile['plasma/1d/ni'][:,0]   
#        ni = np.trapz(ni, rho)
        Zi = self.infile['plasma/znum'][0]
        lnlambda=17
        Ai = self.infile['plasma/anum'][0]
#        ti = self.infile['plasma/1d/ti'][:]*1e-3
#        ti = np.trapz(ti, rho)
        taucoll_i = np.zeros((np.shape(E)[0], np.shape(ni)[0]))
        taucoll_e = np.zeros((np.shape(E)[0], np.shape(ni)[0]))
        self.taucoll_e = np.zeros((np.shape(E)[0], len(np.where(rho<1.)[0])))
        self.taucoll_i = np.zeros((np.shape(E)[0], len(np.where(rho<1.)[0])))
        self.taucoll_e_mean = np.zeros((np.shape(E)[0]))
        self.taucoll_i_mean = np.zeros((np.shape(E)[0]))
        
        for en_ind, en in enumerate(E):
            taucoll_e[en_ind,:] = 1.09e16*en**1.5/(ni*Zi**2*lnlambda)
            taucoll_e[en_ind, rho>=1.] = np.NaN
            self.taucoll_e[en_ind,:] = taucoll_e[en_ind,~np.isnan(taucoll_e[en_ind,:])]
            taucoll_i[en_ind,:] = 6.6e17*Ai**0.5*en**1.5/(ni*Zi**4*lnlambda)
            taucoll_i[en_ind, rho>=1.] = np.NaN
            self.taucoll_i[en_ind,:] = taucoll_i[en_ind,~np.isnan(taucoll_i[en_ind,:])]
            self.taucoll_e_mean[en_ind] = np.trapz(self.taucoll_e[en_ind,:], rho[rho<1.])
            self.taucoll_i_mean[en_ind] = np.trapz(self.taucoll_i[en_ind,:], rho[rho<1.]) 

class TCV_orbits(dat_particles):
    def __init__(self, infile_n, fname_surf=''):
        self.device = 'TCV'
        self.id = infile_n[-12:-4]
        self.a  = 0.25
        self.R0 = 0.88
        dat_particles.__init__(self, infile_n, fname_surf=fname_surf)	


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
