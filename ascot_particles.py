from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, colors
import collections
import math
import h5py
#import ascot_distributions

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

my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)

class particles:
    """
    SUPERCLASS
    
    """
    def __init__(self):
        
        # Initialisation
		self.npart   = 0
		self.nfields = 0
		self.field   = []
		self.unit    = []
		self.data_i  = collections.defaultdict(list)
		self.data_e  = collections.defaultdict(list)
		self.R_w = []
		self.z_w = []
		self._RZsurf()
		self._readwall()

    def _calc_binarea(self,x,y):
        """
        hidden method to calculate the XY area of the bin for XY ionisation plot
        """
        Dx = x[1]-x[0]
        Dy = y[1]-y[0]
        area = Dx*Dy
        return area

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
            
            
    def TCV_calc_shinethr(self):
         """
         Method to calculate the shine-through with the weights
         """

         if self.endgroup != 'shinethr':
             print("WARNING! Check input file, which is ", self.fname)

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

    def plot_RZpart(self, ax=0):
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
        _plot_2d(x, y, xlab=xlab, ylab=ylab, Id=self.id, title='RZ ionization',\
                 wallrz=wallrz, surf=surf, ax=ax)
          
    def plot_XY(self, ax=0):
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
            phi = self.data_i['phiprt']

        if np.mean(R)==999. or np.mean(R)==-999.:
            R=self.data_i['R']
            z=self.data_i['z']
            phi = self.data_i['phi']*math.pi/180.
        x=np.zeros(self.npart)
        y=np.zeros(self.npart)
        for i, el in enumerate(R):
            x[i]=el*math.cos(phi[i])
            y[i]=el*math.sin(phi[i])
        xlab = 'X [m]'
        ylab = 'Y [m]' 
        R0 = 0 
        if self.R0:
            R0 = self.R0
            
        wallxy= [self.R_w, self.z_w]
        _plot_2d(x, y, xlab=xlab, ylab=ylab, Id=self.id, title='XY Ionization',\
                 wallxy=wallxy, R0=R0, ax=ax)

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


    def losses_prediction(self):
        """
        Plots the initial values of the particles depending on their endcondition
        """
        endstate=self.data_e['endcond']
        R   = self.data_i['R']
        phi = self.data_i['phi']
        plt.figure()
        plt.scatter(R*np.cos(phi), R*np.sin(phi))
        ind_th   = endstate == 4
        ind_wall = endstate == 3
        ind_emin = endstate == 2
        plt.figure()
        plt.scatter(self.data_i['pitch'][ind_wall], self.data_i['energy'][ind_wall],c='blue', label='wall', marker='x')
        plt.scatter(self.data_i['pitch'][ind_th],   self.data_i['energy'][ind_th], c='red',  label='th', marker='x')
        plt.scatter(self.data_i['pitch'][ind_emin], self.data_i['energy'][ind_emin],c='green', label='emin', marker='x')
        plt.legend(loc='best')
        
        plt.figure()
        plt.scatter(R[ind_wall]*np.cos(phi[ind_wall]), R[ind_wall]*np.sin(phi[ind_wall]),c='blue', label='wall', marker='x')
        plt.scatter(R[ind_th]*np.cos(phi[ind_th]), R[ind_th]*np.sin(phi[ind_th]),c='red', label='th', marker='x')
        plt.scatter(R[ind_emin]*np.cos(phi[ind_emin]), R[ind_emin]*np.sin(phi[ind_emin]),c='green', label='emin', marker='x')
        
        plt.legend(loc='best')


    def _readwall(self):
        """
        Hidden method to read the wall
        """
        in_w_fname = 'input.wall_2d'
        try:
            wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=1)
        except:
            if self.id[0:2]=='00':
                in_w_fname = '/home/vallar/ASCOT/runs/JT60SA/002/input.wall_2d'
            else:
                in_w_fname = '/home/vallar/ASCOT/runs/TCV/57850/input.wall_2d'
            wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=1)
            
        self.R_w = wall[0,:]
        self.z_w = wall[1,:]
        self.R_w = np.array(self.R_w)
        self.z_w = np.array(self.z_w)

    def _RZsurf(self):
        """
        Reads the position of RZ surfaces from ascot file
        """       
        
        if self.id[0:3]=='003':
            strt=str(self.id[0:3])+'000'
        elif self.id[0:3]=='005':
            strt=str(self.id[0:3])+'053'
        elif self.id[0:5]=='57850' and self.id[5:7]=='08':
            strt =str(self.id[0:5])+'080'
        elif self.id[0:5]=='57850' and self.id[5:7]=='14':
            strt =str(self.id[0:5])+'140'            
        else:
            print('Impossible to load RZ equiflux surfaces')
            print('File ID:'+ self.id+' '+self.id[5:6])
            self.RZsurf = 0
            return 

        f = h5py.File('ascot_'+strt+'.h5')
        print("READING SURF FROM ascot_"+strt+".h5")
        self.RZsurf = f['bfield/2d/psi'].value
        self.Rsurf = f['bfield/r']
        self.zsurf = f['bfield/z']
        edge = f['boozer/psiSepa'][:]; axis=-1*f['boozer/psiAxis'][:]
        self.R0 = f['misc/geomCentr_rz'][0]
        #if self.id[0:3]=='005':
        #    print("Multiply times -1")
        #    self.RZsurf=-1.*self.RZsurf
        self.RZsurf = (self.RZsurf-axis)/(edge-axis)
        self.RZsurf = np.sqrt(self.RZsurf)            
        
    def _plot_RZsurf(self, ax):
        try:
            self.RZsurf.mean()
        except:
            self._RZsurf()
            
        CS = ax.contour(self.Rsurf, self.zsurf, self.RZsurf, [0.2, 0.4, 0.6, 0.8, 1.0], colors='k')
        plt.clabel(CS, inline=1, fontsize=10) 

    def endcondition(self):
        """
        Computes the endcondition of the particles and stores a dict with the amount of MC particles with that final condition
        """
        errendcond  = {'aborted':-2, 'rejected':-1}
        counter = 0
        for key in errendcond:
            ind = np.where(self.data_e['endcond']==errendcond[key])[0]
            if len(ind)!=0:
                counter += len(ind)
                print("STRANGE END CONDITION! {} ({:d}) : {:d} particles".format(key, errendcond[key], len(ind)))

        physcond = {'none':0,'tmax':1,'emin':2,'wall':3,'thermalisation':4}
        for key in physcond:
            ind = np.where(self.data_e['endcond']== physcond[key])[0]
            counter += len(ind)
            print("{} ({:2d}) : {:4d} particles".format(key, physcond[key], len(ind)))

        infcond = {'cputmax':17, 'outsidemaxrho':10, 'outsidewall':16}
        for key in infcond:
            ind = np.where(self.data_e['endcond']== infcond[key])[0]
            counter += len(ind)
            print("{} ({:2d}) : {:4d} particles".format(key, infcond[key], len(ind)))    

        print("Total particles counted :", counter/len(self.data_e['endcond'])*100., ' %')

    def particle_status(self):
        """
        With this function you can plot the histogram with thermalised and lost particles
        """
        # text font for plot
        plt.rc('font', family='serif', serif='Palatino')
        #plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=========================

        width=0.8

        fig=plt.figure(figsize=(12,5))
        ax=fig.add_subplot(111)
        ax.bar(np.arange(14)+0.5,self.endstatus[:,2], width, color='g', label='Th.')
        if n_wall!=0:      
            ax.bar(np.arange(14)+0.5,self.endstatus[:,1], width, bottom=self.endstatus[:,2], color='r', label='Wall')
            frac_wall = self.endstatus[:,1]/self.npart_arr[:]
            frac_tot_wall = np.sum(self.endstatus[:,1])/np.sum(self.npart_arr)
        ax.bar(np.arange(14)+0.5,self.endstatus[:,3], width, bottom=self.endstatus[:,2]+self.endstatus[:,1], color='y', label='$E_{min}$')        
        if n_rho!=0:        
            ax.bar(np.arange(14)+0.5,self.endstatus[:,4], width, bottom=self.endstatus[:,3]+self.endstatus[:,2]+self.endstatus[:,1], color='b', label=r'$\rho_{MAX}$')        
            frac_rho = self.endstatus[:,4]/self.npart_arr[:]
            frac_tot_rho = np.sum(self.endstatus[:,4])/np.sum(self.npart_arr)

        # Create your ticker object with M ticks
        M = 4
        yticks = ticker.MaxNLocator(M)
        # Set the yaxis major locator using your ticker object. You can also choose the minor
        # tick positions with set_minor_locator.
        ax.yaxis.set_major_locator(yticks)

        #SET TICKS LABELS ON X AXIS
        ax.xaxis.set_ticks(['NBI'])
        ax.xaxis.set_ticklabels(xlab)

        
        plt.ylabel('Particles')
        plt.xlabel('Beam index')
        plt.legend(loc='best', prop={'size':20})

        print("Fraction to the wall:",frac_wall)
        print("Fraction outside the rho:", frac_rho)
        print("Total fraction to wall:", frac_tot_wall)
        print("Total fraction outside rho:", frac_tot_rho)

    def plot_initial_wall_pos(self, ax=0):
        """
        Plot with final position of the particles colliding with wall and energy 
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        r = self.data_i['R'][ind]
        z = self.data_i['z'][ind]
        R0=self.infile['misc/geomCentr_rz'][0]
        theta = np.arctan2(z,r-R0)
        energy = self.data_e['energy'][ind]*1e-3
        pitch = self.data_i['pitch'][ind]
        ind = np.where(np.abs(theta)<0.4)
        r = r[ind]
        z = z[ind]
        theta=theta[ind]
        phi = self.data_i['phi'][ind]
        wallrz= [self.R_w, self.z_w]
        energy = energy[ind]
        pitch = pitch[ind]

        
        _plot_2d(r,z, 'R [m]', 'z [m]', scatter=energy, Id=self.id, wallrz=wallrz, surf=[self.Rsurf, self.zsurf, self.RZsurf], ax=ax)
        ax=0
        ax=plt.gca()
        ax.axvline(x=np.max(self.R_w-0.1))
        plt.figure()
        plt.hist(pitch, bins=20)        
        _plot_2d(phi, theta, xlab=r'$\phi$',ylab=r'$\theta$', Id = self.id, hist=1)

    def plot_wall_losses(self, ax=0):
        """
        Plot with final position of the particles colliding with wall and energy 
        """
        ind = np.where(self.data_e['endcond']== 3)[0] #wall
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        R0=self.infile['misc/geomCentr_rz'][0]
        theta = np.arctan2(z,r-R0)
        phi = self.data_e['phi'][ind]
        wallrz= [self.R_w, self.z_w]
        energy = self.data_e['energy'][ind]*1e-3

        _plot_2d(r,z, 'R [m]', 'z [m]', scatter=energy, Id=self.id, wallrz=wallrz, surf=[self.Rsurf, self.zsurf, self.RZsurf], axin=ax)
        
        #_plot_2d(phi, theta, xlab=r'$\phi$',ylab=r'$\theta$', Id = self.id, hist=1)

    def maxrho_histo(self):
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
        pitch = self.data_e['pitch'][ind]
        vr = self.data_e['vR'][ind]
        vz = self.data_e['vz'][ind]
        vphi = self.data_e['vphi'][ind]
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        R0=self.infile['misc/geomCentr_rz'][0]
        theta = np.arctan2(z,r-R0)
        phi = self.data_e['phi'][ind]
        x = r*np.cos(phi); y=r*np.sin(phi)        
      
        
        
        #plt.close('all')
#        plt.figure(); plt.hist(pitch, bins=20); plt.xlabel('Pitch'); plt.ylabel('Number of particles')
#        plt.figure(); plt.hist(energy, bins=30); plt.xlabel('Energy'); plt.ylabel('Number of particles')
##        plt.figure(); plt.hist(vr*1e-3, bins=20); plt.xlabel(r'$v_r$ [km/s]'); plt.ylabel('Number of particles')
##        plt.figure(); plt.hist(vz*1e-3, bins=20); plt.xlabel(r'$v_z$ [km/s]'); plt.ylabel('Number of particles')
#        plt.figure(); plt.hist(phi, bins=20); plt.xlabel('Phi (toroidal angle)'); plt.ylabel('Number of particles')
#        plt.figure(); plt.hist(theta, bins=20); plt.xlabel('theta (poloidal angle)'); plt.ylabel('Number of particles')
#



        plt.figure(); plt.scatter(x,y);  plt.grid('on'); plt.xlabel(r'x'); plt.ylabel('y')
#       
#        plt.figure(); plt.scatter(pitch, energy*1e-3);
#        plt.scatter(pitchi, energyi*1e-3, color='r');  plt.grid('on'); plt.xlabel(r'Pitch $v_\parallel/v$'); plt.ylabel('Energy [keV]')


        #plt.figure(); plt.scatter(vphi*1e-3, np.sqrt(vr**2+vz**2)*1e-3);  plt.grid('on'); plt.xlabel(r'$v_\parallel$ [km/s]'); plt.ylabel(r'$v_\perp [km/s]$');
        #plt.figure(); plt.scatter(vr*1e-3, vz*1e-3);  plt.grid('on'); plt.xlabel(r'$v_r$ [km/s]'); plt.ylabel(r'$v_z$ [km/s]'); plt.axis('equal')

        #plt.figure(); plt.scatter(r,z, c=pitch); plt.colorbar(); plt.xlabel('R'); plt.ylabel('z')
        
        #ax=plt.axes(); ax.quiver(r,z,np.multiply(vr,1./norm_v), np.multiply(vz, 1./norm_v));
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
        endcond = self.data_e['endcond']        
        ind = np.where(np.logical_or(endcond == 1, endcond == 2))  #Tmax, Emin, cputmax
        p_res = np.dot(self.data_e['weight'][ind],self.data_e['energy'][ind])*1.602e-19
        self.pcoup = p_ini-p_end+p_res

class dat_particles(particles):
    """
    Class (inherited from particles) handling dat files (e.g. input.particles, test_ascot.orbits.dat)
    """
    def __init__(self, infile_n):
        """
        READ ASCII FILE
        """
        particles.__init__(self)
        self.fname = infile_n
        in_f  = open(self.fname)
        lines = in_f.readlines()[3:]
        in_f.close()
        self.id = self.fname[-10:-4]
        self.infile = h5py.File('ascot_'+self.id+'.h5')
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
            self.npart=50000
            #self._compute_npart(lines[-2])

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
        for i in range(self.npart):
            for key in self.data_i:
                self.data_i[key] = np.append(self.data_i[key], \
                            self.partdict[i][key][0])
#        for key in self.field:
#            if key=='id':
#                self.data_i[key] = np.array(self.data_i[key], dtype=int)
#            else:
#                self.data_i[key] = np.array(self.data_i[key], dtype=float)

        #CONVERT PHI FROM DEG TO RAD
        #self.origins = np.array(list(set(self.data_i['origin'])))
        print("FIELDS IN FILE ",self.fname," :")
        print(self.field)
        self.banana_orbit()

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
        self._banana_dim()
        
    def _banana_dim(self):
        """
        computes the fraction of trapped particles with the following formula:
         ft = sqrt(epsilon)*(1.46-0.46*epsilon)
         
        and the width of the banana with the following formula
         w_b = sqrt(epsilon)*rho_L(larmor radius evaluated using only poloidal field)
        """
        try:
            self.partdict[0]['mass'].mean()
            self.partdict[0]['charge'].mean()
        except:
            print("Impossible to calculate the banana orbits dimension")
            return
        try:
            self.partdict[0]['charge'].mean()
        except:
            for i in self.partdict:
                i['charge']=1
            
        R_torus = 2.96 #in meters
        a = 1.11
        self.epsilon = a/R_torus
        print("Epsilon used is for scenario 5")      
        #Calculating larmor radius using only poloidal field
        self.rhoL_poloidal = np.zeros(self.npart, dtype=float)
        for i in range(self.npart):
            if self.trappind[i]==0:
                continue
            m,q,E = self.partdict[i]['mass'], self.partdict[i]['charge'], \
                    self.partdict[i]['energy']
            q = q*1.602e-19; E = E*1.602e-19
            m = m*1.66e-27
            m = m[0]; q=q[0]; E=E[0]
            v = (2.*E/m)**0.5
            Bpol = np.sqrt(self.partdict[i]['BR']**2+self.partdict[i]['Bz']**2)
            Bpol = Bpol[0]            
            #print m,v,q,Bpol
            self.rhoL_poloidal[i] = (m*v)/(q*Bpol)
        
        self.w_b_unfilt = self.rhoL_poloidal*self.epsilon**0.5
        self.w_b = self.w_b_unfilt[self.w_b_unfilt<1.]
        
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
        plt.rc('font', size=20)

        ylim=[-1., -0.5]
        
        self._RZsurf()
        ind_t = np.where(self.trappind==1)[0]
        ind_p = np.where(self.trappind==0)[0]     
        CL = ['','']
        r = self.data_i['R']; z = self.data_i['z']
        p = self.data_i['pitch']
        rho = self.data_i['rho']
        qpart = np.zeros(self.npart)
#        q = self.infile['/boozer/qprof'][:]
#        rho_q = self.infile['/boozer/psi'][:]**0.5
#        for i in range(self.npart):
#            ind = np.argmin(rho[i]-rho_q>0)
#            qpart[i] = q[ind]
            
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
        
        self._plot_RZsurf(axrz)
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

        self._plot_RZsurf(axrz2)
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
        
    def plot_trapped_energy_PNBs(self):
        """
        plots trapped particles as function of energy
        """
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
        f.suptitle(self.id)
        ax = f.add_subplot(111)
        f.text(0.01, 0.95, str(float(self.ntrapp)/self.npart))
        x = [85000./3.,85000./2.,85000.]
        x = [0.33-width/2., 0.66-width/2., 1.-width/2.]
        y = [float(num_t_thir)/num_thir, float(num_t_half)/num_half, float(num_t_full)/num_full]
        y1 = [float(num_thir)/num_tot, float(num_half)/num_tot, float(num_full)/num_tot]
        y2 = [float(num_t_thir)/num_tot, float(num_t_half)/num_tot, float(num_t_full)/num_tot]

        ax.bar(x, y1, width, color='b', label='Passing')
        ax.bar(x, y2, width, color='r', label='Trapped')
        ax.set_ylim([0, 1.4])
        ax.set_xticks([0.34, 0.66, 1.])
        ax.set_xticklabels([r'$E_{inj}$/3', r'$E_{inj}$/2', r'$E_{inj}$'])
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.])
        
        for i in range(3):
            ax.text(x[i]+0.05, y1[i]+0.05, '{:.1f} %'.format(y[i]*100.), color='r')

        ax.legend(loc='best')
        ax.set_ylabel(r'Fraction to total number')
        ax.yaxis.grid('on')

        
    def plot_trapped_energy_NNBs(self):
        """
        plots trapped particles as function of energy
        """
        e = self.data_i['energy']; 
        width=0.1
        num_t_full = len(np.where(self.trappind==1)[0])
        num_full = len(e)

        f = plt.figure()
        f.suptitle(self.id)
        ax = f.add_subplot(111)
        x = [500e3]
        x = [0.5]
        y = [float(num_t_full)/num_full]

        ax.bar(x, [1.], width, color='b', label='Passing')
        ax.bar(x, y, width, color='r', label='Trapped')
        ax.set_ylim([0, 1.4])
        ax.set_xlim([0.45, 0.65])
        ax.legend(loc='best')
        ax.yaxis.grid('on')
        ax.set_xticks([0.55])
        ax.set_xticklabels([r'$E_{inj}$'])
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.])
        ax.set_ylabel(r'Fraction to total number')
        ax.text(x[0], 1.05, '{:.1f} %'.format(y[0]*100.), color='r')
        
    def plot_trajectory(self):
        """
        Plots trajectory of the particles
        """
        try:
            self.trappind.mean()
        except:
            self.banana_orbit()
            
        ind_t = np.where(self.trappind==1)[0]
        ind_p = np.where(self.trappind==0)[0]      
        ind = np.linspace(0, len(ind_t)-1, len(ind_t)-1)
        f = plt.figure()
        f.suptitle(self.fname)
        axtrajRZ = f.add_subplot(121)        
        for i in ind:
            axtrajRZ.plot(self.data_i['R'][ind_t[i]], self.data_i['z'][ind_t[i]], 'k*')
            axtrajRZ.plot(self.partdict[ind_t[i]]['R'], self.partdict[ind_t[i]]['z'], 'kx') 
            #axtrajRZ.plot(self.data_i['R'][ind_p[i]], self.data_i['z'][ind_p[i]], 'r*')
            #axtrajRZ.plot(self.partdict[ind_p[i]]['R'], self.partdict[ind_p[i]]['z'], 'r-') 
        #self._plot_RZsurf(axtrajRZ)
        axtrajRZ.set_xlabel(r'R (m)')
        axtrajRZ.set_ylabel(r'z (m)')
        
        axtrajxy = f.add_subplot(122)        
#         for i in ind:
#             ll = ind_t[i]
#             xi,yi = self.data_i['R'][ll]*np.cos(self.data_i['phi'][ll]), \
#                     self.data_i['R'][ll]*np.sin(self.data_i['phi'][ll])
#             x,y = self.partdict[ll]['R']*np.cos(self.partdict[ll]['phi']), \
#                     self.partdict[ll]['R']*np.sin(self.partdict[ll]['phi'])  
                    
#             axtrajxy.plot(xi,yi, 'k*')
#             axtrajxy.plot(x,y, 'k-') 
#             ll = ind_p[i]
#             xi,yi = self.data_i['R'][ll]*np.cos(self.data_i['phi'][ll]), \
#                     self.data_i['R'][ll]*np.sin(self.data_i['phi'][ll])
#             x,y = self.partdict[ll]['R']*np.cos(self.partdict[ll]['phi']), \
#                     self.partdict[ll]['R']*np.sin(self.partdict[ll]['phi'])             
#             axtrajxy.plot(xi,yi, 'r*')
#             axtrajxy.plot(x,y, 'r-') 

        axtrajxy.set_xlabel(r'x (m)')
        axtrajxy.set_ylabel(r'y (m)')        

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
       
class h5_particles(particles):
	"""
	Class (inherited from particles) handling h5 files (e.g. bbnbi.h5, ascot.h5)
	"""
	def __init__(self, fname):
		"""
		Initialising
		"""
		self.fname = fname
		indd = self.fname[-9:-3]
		self.id = indd
		if indd[0:2] != '00':
			self.id = self.fname[-11:-3]
			
		particles.__init__(self)
		dataf = h5py.File(self.fname)
		self.infile=dataf
		try:
			self.endgroup = 'shinethr' #this is for bbnbi
			self.field = dataf[self.endgroup].keys()
		except:
			self.endgroup = 'endstate' #this is for ascot
			self.field = dataf[self.endgroup].keys()
			
		self.nfields=len(self.field)      
		for key in self.field:
			if key=='id':
				self.data_i[key] = np.array(dataf['inistate/'+key], dtype=int)
				self.data_e[key] = np.array(dataf[self.endgroup+'/'+key], dtype=int)
			else:
				self.data_i[key] = np.array(dataf["inistate/"+key], dtype=float)
				self.data_e[key] = np.array(dataf[self.endgroup+"/"+key], dtype=float)

		if self.endgroup == 'endstate':
			self._h5origins = dataf['species/testParticle/origin'].value
		elif self.endgroup == 'shinethr':
			self.origins = np.sort(np.array(list(set(self.data_i['origin']))))

		# evaluate number of particles
		self.npart = np.max(self.data_i['id'])

		#CONVERT PHI FROM DEG TO RAD
		self.data_i['phiprt']=self.data_i['phiprt']*math.pi/180.
		self.data_e['phiprt']=self.data_e['phiprt']*math.pi/180.

		print("FIELDS IN FILE ",self.fname,", section inistate and ",self.endgroup," :")
		print(list(self.field))

		# Swapping R and Rprt, since they are opposite to matlab
		self.data_e['R']   = self.data_e['Rprt']
		self.data_e['z']   = self.data_e['zprt']
		self.data_e['phi'] = self.data_e['phiprt']

	def calc_shinethr(self):
		"""
		Method to calculate the shine-through with the weights
		"""
	#        try:
	#            self.originpower.mean()
	#        except:
	#            self._calc_originbeam()
			
			
		if self.endgroup != 'shinethr':
			print("WARNING! Check input file, which is ", self.fname)

		self._calc_weight()
		self.shinethr=np.zeros((26), dtype=float)
		self.shinethr_abs=np.zeros((26), dtype=float)
		self.id2beamnum = \
						{\
						'45':1 ,  '46':1,    '47':2,    '48':2,   \
						'133':3,  '134':3,   '135':4,   '136':4,  \
						'221':5,  '222':5,   '223':6,   '224':6,  \
						'309':7,  '310':7,   '311':8,   '312':8,  \
						'3637':9, '3638':9,  '3639':10, '3640':10,\
						'5253':13,'5254':13, '5255':14, '5256':14,\
						'3031':99,'3032':101 \
						}
		for ind_i,i in enumerate(self.id2beamnum):
			ind = self.data_e['origin'][:]==int(i)
			if len(self.data_e['origin'][ind])==0:
				w=0
				e=0
				power = 1
			else:  
				power=1.e6
				if int(i) in [3031, 3032]:
					power = 5.e6
					
				e = self.data_e['energy'][ind][0]
				e = e*1.612e-19
				w = np.sum(self.data_e['weight'][ind_i])
				#wtot = self.w[ind_i]
			self.shinethr[ind_i]=float(e)*w/power
			self.shinethr_abs[ind_i] = float(e)*w


def _plot_2d(x, y, xlab, ylab, Id='', title='', wallxy=0, wallrz=0, surf=0, R0=0, axin=0, scatter=0, hist=0,  **kwargs):
	"""
	Hidden method to plot the 2D histogram
	wall: set to 1 if wall needed to plot (i.e. RZ function)
	"""
	#==============================================
	# SET TEXT FONT AND SIZE
	#==============================================
	#plt.rc('font', family='serif', serif='Palatino')
	#plt.rc('text', usetex=True)
	plt.rc('xtick', labelsize=20)
	plt.rc('ytick', labelsize=20)
	plt.rc('axes', labelsize=20)
	#==============================================
	
	flag_dict = kwargs
	figsize=[8,8]; flag_label=1

	if axin==0:
            if wallrz != 0 :
                figsize=[6,7]
            # Defining figure and ax
            fig = plt.figure(figsize=figsize)
            fig.text(0.01, 0.01, Id)
            ax  = fig.add_subplot(111)
        else:
            ax=axin
            flag_label=0
            fig = plt.gcf()

        #doing the actual plot
        if len(fig.axes)==1:
            or_cb = 'vertical'
        else:
            or_cb = 'horizontal'
        
        if np.mean(scatter)!=0:
            pp=ax.scatter(x, y, 100, c=scatter)
            plt.colorbar(pp, ax=ax, orientation = or_cb)
        elif hist != 0:
            ax.hist2d(x, y, bins=30, cmap=my_cmap)
        else:
            hb = ax.hist2d(x, y, bins=100, cmap=my_cmap)
            fig.colorbar(hb[3], ax=ax, orientation=or_cb)

        
	#Checks for wall and plots it	
	if wallrz != 0:
            ax.plot(wallrz[0], wallrz[1], 'k', linewidth=3)
            ax.axis('equal')
	elif wallxy != 0:
            rmin = np.min(wallxy[0])
            rmax = np.max(wallxy[0])
            circlemin = plt.Circle((0,0), rmin, color='k', fill=False, linewidth=3)
            circlemax = plt.Circle((0,0), rmax, color='k', fill=False, linewidth=3)
            ax.add_artist(circlemin); ax.add_artist(circlemax)
            ax.axis('equal')
	# Checks for magnetic axis plotting
	if R0!=0:
            circle1 = plt.Circle((0, 0), R0, color='r', fill=False, linestyle='--')      
            ax.add_artist(circle1)
	#Checks for magnetic surfaces and plots them
	if surf!= 0 and axin==0:
            try:
                llines = [0.2, 0.4, 0.6, 0.8, 1.0]
                CS = ax.contour(surf[0], surf[1], surf[2], llines, colors='k')
                plt.clabel(CS, inline=1, fontsize=10) 
            except:
                print("Impossible to plot RZ surfaces")     
		
	if flag_label==1 :
            #==================================
            # SET TICK LOCATION
            #==================================

            # Create your ticker object with M ticks
            M = 4
            yticks = ticker.MaxNLocator(M)
            xticks = ticker.MaxNLocator(M)
            # Set the yaxis major locator using your ticker object. You can also choose the minor
            # tick positions with set_minor_locator.
            ax.yaxis.set_major_locator(yticks)
            #ax.yaxis.set_minor_locator(yticks_m)
            ax.xaxis.set_major_locator(xticks)
            #==================================
	
            ax.set_title(title)
            ax.set_xlabel(xlab);ax.set_ylabel(ylab)	
            ax.grid('on')

        fig.tight_layout()
        plt.show()
	if 'fname' in flag_dict:
		plt.savefig(flag_dict['fname'], bbox_inches='tight')
