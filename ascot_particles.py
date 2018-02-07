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
            print "No weights available"
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
            print "WARNING! Check input file, which is ", self.fname

        self._calc_weight()
        self.shinethr=np.zeros((1), dtype=float)
        self.shinethr_abs=np.zeros((1), dtype=float)
        power = 1e6

        e = self.data_e['energy']
        e = e*1.612e-19
        w = self.data_e['weight']
        self.shinethr_abs= np.sum(e*w)
        self.shinethr=np.sum(e*w)/power



        print "TCV Shine-through:", "{0:5f} %".format(self.shinethr*100),\
                  "||  {0:3f} W".format(self.shinethr_abs)


            
    def plot_RZpart(self):
        """
        Method to plot R vs z of the ionised particles, useful mostly with bbnbi
        """
        try:
            x=self.data_i['Rprt']
            y=self.data_i['zprt']
        except:
            x=self.data_i['R']
            y=self.data_i['z']
        if np.mean(x)==999.:
            x=self.data_i['R']
            y=self.data_i['z']
            
        f=plt.figure()
        f.suptitle(self.fname)
        ax = f.add_subplot(111)
        self._plot_RZsurf(ax)
        #hb = ax.hexbin(x, y, gridsize=100, cmap=my_cmap)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap)
        f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('R')
        ax.set_ylabel('z')

        if len(self.R_w)==0:
            self._readwall()
        ax.plot(self.R_w , self.z_w, 'k', linewidth=3)

        ax.axis('equal')
        ax.axis([min(self.R_w), max(self.R_w), min(self.z_w), max(self.z_w)])
        #ax.set_xrange([min(self.R_w), max(self.R_w)])
        #ax.set_yrange([min(self.z_w), max(self.z_w)])                   
        plt.show()


    def plot_XYpart_singlebeam(self):
        """
        Method to plot XY of ionisation, without difference between the beams
        """

        x=np.zeros(self.npart)
        y=np.zeros(self.npart)
        
        for i, el in enumerate(self.data_i['Rprt']):
            x[i]=el*math.cos(self.data_i['phiprt'][i])
            y[i]=el*math.sin(self.data_i['phiprt'][i])
            
        f=plt.figure()
        ax=f.add_subplot(111)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap, normed=True)
        f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        theta=np.arange(0,6.3,0.02*6.28)
        if len(self.R_w)==0:
            self._readwall()
        ax.plot(np.min(self.R_w)*np.cos(theta) , np.min(self.R_w)*np.sin(theta), 'k', linewidth=3)
        ax.plot(np.max(self.R_w)*np.cos(theta) , np.max(self.R_w)*np.sin(theta), 'k', linewidth=3)
        plt.show()

    def plot_ioniz_beams(self):
        """
        Method to plot XY of ionisation FOR EACH BEAM
        """
        try:
            self.beamorigindict.keys()
        except:
            self._calc_originbeam()

            
        #x_range=[-5,5]
        axxy = plt.subplot2grid((1,2),(0,0))
        axrz = plt.subplot2grid((1,2),(0,1))
        ind = np.zeros((2,self.npart), dtype=bool)
        
        for i,el in enumerate(self.beamorigindict.keys()):
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
        
##         axxy.('X')
##         axxy.ylabel('Y')
##         axrz.xlabel('R')
##         axrz.ylabel('z')
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
            in_w_fname = '/home/vallar/ASCOT/runs/JT60SA/002/input.wall_2d'
            wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=1)
        self.R_w = wall[0,:]
        self.z_w = wall[1,:]
        self.R_w = np.array(self.R_w)
        self.z_w = np.array(self.z_w)

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

        print "Total particles counted :", counter/len(self.data_e['endcond'])*100., ' %'

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

        print "Fraction to the wall:",frac_wall
        print "Fraction outside the rho:", frac_rho
        print "Total fraction to wall:", frac_tot_wall
        print "Total fraction outside rho:", frac_tot_rho

    def maxrho_histo(self):
        """
        Histogram with final theta position, pitch, energy and 2D plot of the particle velocity
        """
        try:
            self.R_w.mean()
        except:
            self._readwall()
        #ind = np.where(self.data_e['endcond']== 10)[0] #rho=1
        #ind = np.where(self.data_e['endcond']== 3)[0] #wall
        ind = np.arange(len(self.data_e['endcond'])) #all particles
        
        pitchi = self.data_i['pitch'][ind]
        energyi = self.data_i['energy'][ind]
        pitch = self.data_e['pitch'][ind]
        energy = self.data_e['energy'][ind]
        vr = self.data_e['vR'][ind]
        vz = self.data_e['vz'][ind]
        vphi = self.data_e['vphi'][ind]
        r = self.data_e['R'][ind]
        z = self.data_e['z'][ind]
        theta = np.arctan2(z,r-3.)*180./math.pi
        phi = self.data_e['phi'][ind]
        x = r*np.cos(phi); y=r*np.sin(phi)        
            
        
        plt.close('all')
        plt.figure(); plt.hist(pitch, bins=20); plt.xlabel('Pitch'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(energy, bins=30); plt.xlabel('Energy'); plt.ylabel('Number of particles')
#        plt.figure(); plt.hist(vr*1e-3, bins=20); plt.xlabel(r'$v_r$ [km/s]'); plt.ylabel('Number of particles')
#        plt.figure(); plt.hist(vz*1e-3, bins=20); plt.xlabel(r'$v_z$ [km/s]'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(phi, bins=20); plt.xlabel('Phi (toroidal angle)'); plt.ylabel('Number of particles')
        plt.figure(); plt.hist(theta, bins=20); plt.xlabel('theta (poloidal angle)'); plt.ylabel('Number of particles')

        plt.figure(); plt.plot(self.R_w, self.z_w, 'k', lw=2.3); plt.scatter(r,z, 100, c=energy); plt.title('Energy');plt.colorbar(); plt.xlabel('R'); plt.ylabel('z')
        plt.figure(); plt.scatter(x,y);  plt.grid('on'); plt.xlabel(r'x'); plt.ylabel('y')
       
        plt.figure(); plt.scatter(pitch, energy*1e-3);
        plt.scatter(pitchi, energyi*1e-3, color='r');  plt.grid('on'); plt.xlabel(r'Pitch $v_\parallel/v$'); plt.ylabel('Energy [keV]')

        #plt.figure(); plt.scatter(vphi*1e-3, np.sqrt(vr**2+vz**2)*1e-3);  plt.grid('on'); plt.xlabel(r'$v_\parallel$ [km/s]'); plt.ylabel(r'$v_\perp [km/s]$');
        #plt.figure(); plt.scatter(vr*1e-3, vz*1e-3);  plt.grid('on'); plt.xlabel(r'$v_r$ [km/s]'); plt.ylabel(r'$v_z$ [km/s]'); plt.axis('equal')

        #plt.figure(); plt.scatter(r,z, c=pitch); plt.colorbar(); plt.xlabel('R'); plt.ylabel('z')
        
        #ax=plt.axes(); ax.quiver(r,z,np.multiply(vr,1./norm_v), np.multiply(vz, 1./norm_v));
        
        
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

    def _RZsurf(self):
        """
        Reads the position of RZ surfaces from ascot file
        now the edge is set to the value for scenario 5 from JT60SA
        """       
        f = h5py.File('ascot_'+self.id+'.h5')
        self.RZsurf = f['bfield/2d/psi'].value
        self.Rsurf = f['bfield/r']
        self.zsurf = f['bfield/z']
        edge = np.abs(f['boozer/psiSepa'][:]); axis=np.abs(f['boozer/psiAxis'][:])
        self.RZsurf = (self.RZsurf - axis )/(edge-axis)
        self.RZsurf = np.sqrt(self.RZsurf)            

    def _plot_RZsurf(self, ax):
        try:
            self.RZsurf.mean()
        except:
            self._RZsurf()
            
        CS = ax.contour(self.Rsurf, self.zsurf, self.RZsurf, [0.2, 0.4, 0.6, 0.8, 1.0], colors='k')
        plt.clabel(CS, inline=1, fontsize=10)        
        

class dat_particles(particles):
    """
    Class (inherited from particles) handling dat files (e.g. input.particles, test_ascot.orbits.dat)
    """
    def __init__(self, fname):
        """
        READ ASCII FILE
        """
        particles.__init__(self)
        self.fname = fname
        in_f  = open(self.fname)
        lines = in_f.readlines()[3:]
        in_f.close()
        self.id = self.fname[-10:-4]
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
            self._compute_npart(lines[-2])

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
                print "Doing particle ", part_id
            for i, el in enumerate(row):
                if self.field[i]=='phi':
                    self.partdict[part_id][self.field[i]] = \
                    np.append(self.partdict[part_id][self.field[i]], el*math.pi/180.)                   
                else:
                    self.partdict[part_id][self.field[i]] = \
                    np.append(self.partdict[part_id][self.field[i]], el)                

        self.data_i = dict.copy(tmpdict)
        for i in range(self.npart):
            for key in self.data_i.keys():
                self.data_i[key] = np.append(self.data_i[key], \
                            self.partdict[i][key][0])
#        for key in self.field:
#            if key=='id':
#                self.data_i[key] = np.array(self.data_i[key], dtype=int)
#            else:
#                self.data_i[key] = np.array(self.data_i[key], dtype=float)

        #CONVERT PHI FROM DEG TO RAD
        #self.origins = np.array(list(set(self.data_i['origin'])))
        print "FIELDS IN FILE ",self.fname," :"
        print self.field
        

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
            print "Impossible to calculate the banana orbits dimension"
            return
            
        R_torus = 2.96 #in meters
        a = 1.11
        self.epsilon = a/R_torus
        print "Epsilon used is for scenario 5"        
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
        self._RZsurf()
        ind_t = np.where(self.trappind==1)[0]
        ind_p = np.where(self.trappind==0)[0]     
        CL = ['','']
        r = self.data_i['R']; z = self.data_i['z']
        p = self.data_i['pitch']
        rho = self.data_i['rho']
        nbins=50; nsurf=20
        
        hist_p, xedges_p, yedges_p = np.histogram2d(r[ind_p], z[ind_p], bins=nbins, normed=True)
        hist_t, xedges_t, yedges_t = np.histogram2d(r[ind_t], z[ind_t], bins=nbins, normed=True)
        vmaxt = max(np.max(hist_t), np.max(hist_p))
        ind_cb=1
        if np.max(hist_p)> np.max(hist_t):
            ind_cb=0
            
        f=plt.figure()
        f.suptitle(self.id)
        
        axrz = f.add_subplot(221)
        x = np.linspace(np.min(r[ind_p]), np.max(r[ind_p]), num=nbins)
        y = np.linspace(np.min(z[ind_p]), np.max(z[ind_p]), num=nbins)
        CL[0]=axrz.contourf(x,y,hist_p.T, nsurf, cmap=my_cmap, vmin=0, vmax=vmaxt)
        self._plot_RZsurf(axrz)
        axrz.set_title('Passing n(R,Z)'); axrz.set_xlabel('R [m]'); axrz.set_ylabel('Z [m]')        
        axrz.set_xlim([np.min(x), np.max(x)]); axrz.set_ylim([np.min(y), np.max(y)])
        
        axrz2 = f.add_subplot(222)
        x = np.linspace(np.min(r[ind_t]), np.max(r[ind_t]), num=nbins)
        y = np.linspace(np.min(z[ind_t]), np.max(z[ind_t]), num=nbins)
        CL[1]=axrz2.contourf(x,y,hist_t.T, nsurf, cmap=my_cmap, vmin=0, vmax=vmaxt)
        self._plot_RZsurf(axrz2)
        axrz2.set_title('Trapped n(R,Z)'); axrz2.set_xlabel('R [m]'); axrz2.set_ylabel('Z [m]')
        axrz2.set_xlim([np.min(x), np.max(x)]); axrz2.set_ylim([np.min(y), np.max(y)])
        cbar_ax = f.add_axes([0.85, 0.6, 0.03, 0.3])
        f.colorbar(CL[ind_cb], cax=cbar_ax)


        hist_t, yedges_t, xedges_t = np.histogram2d(rho[ind_t], p[ind_t], bins=nbins, normed=True)
        hist_p, yedges_p, xedges_p = np.histogram2d(rho[ind_p], p[ind_p], bins=nbins, normed=True)
        vmaxt = max(np.max(hist_t), np.max(hist_p))   

        f2=plt.figure()
        f2.suptitle(self.id)
        ax2=f2.add_subplot(111)        
        axrp = f.add_subplot(223)
        x = np.linspace(np.min(rho[ind_p]), np.max(rho[ind_p]), num=nbins)
        y = np.linspace(np.min(p[ind_p]), np.max(p[ind_p]), num=nbins)
        CL[0]=axrp.contourf(x,y,hist_p.T, nsurf, cmap=my_cmap, vmin=0, vmax=vmaxt)

        ax2.contour(x,y,hist_p.T,nsurf,colors='k', label='Passing')
        
        axrp.set_title(r'Passing n($\rho$, $\xi$)'); axrp.set_xlabel(r'$\rho$'); axrp.set_ylabel(r'$\xi$')  
        axrp.set_xlim([np.min(x), np.max(x)]); axrp.set_ylim([np.min(y), np.max(y)])
        
        axrp2 = f.add_subplot(224)
        x = np.linspace(np.min(rho[ind_t]), np.max(rho[ind_t]), num=nbins)
        y = np.linspace(np.min(p[ind_t]), np.max(p[ind_t]), num=nbins)
        CL[1] = axrp2.contourf(x,y,hist_t.T, 10, cmap=my_cmap, vmin=0, vmax=vmaxt)
        ax2.contour(x,y,hist_t.T,nsurf,colors='r', label='Trapped')
        axrp2.set_title(r'Trapped n($\rho$, $\xi$)'); axrp2.set_xlabel(r'$\rho$'); axrp2.set_ylabel(r'$\xi$')  
        axrp2.set_xlim([np.min(x), np.max(x)]); axrp2.set_ylim([np.min(y), np.max(y)])        
        ax2.set_xlabel(r'$\rho$'); ax2.set_ylabel(r'$\xi$')
        ax2.legend(loc='upper right')
        
        f.tight_layout()
        f.subplots_adjust(right=0.8)
        cbar_ax = f.add_axes([0.85, 0.1, 0.03, 0.35])
        f.colorbar(CL[ind_cb], cax=cbar_ax)
        plt.grid('on')


    def plot_trapped_energy(self):
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
        
        f = plt.figure()
        f.suptitle(self.id)
        ax = f.add_subplot(111)
        x = [85000./3.,85000./2.,85000.]
        x = [0.33-width/2., 0.66-width/2., 1.-width/2.]
        y = [float(num_t_thir)/num_thir, float(num_t_half)/num_half, float(num_t_full)/num_full]

        ax.bar(x, [1., 1., 1.], width, color='b', label='Passing')
        ax.bar(x, y, width, color='r', label='Trapped')
        ax.set_ylim([0, 1.4])
        ax.legend(loc='best')
        ax.grid('on')
        
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
        ind = np.linspace(0, 100, num=5)
        f = plt.figure()
        f.suptitle(self.fname)
        axtrajRZ = f.add_subplot(121)        
        for i in ind:
            axtrajRZ.plot(self.data_i['R'][ind_t[i]], self.data_i['z'][ind_t[i]], 'k*')
            axtrajRZ.plot(self.partdict[ind_t[i]]['R'], self.partdict[ind_t[i]]['z'], 'k-') 
            axtrajRZ.plot(self.data_i['R'][ind_p[i]], self.data_i['z'][ind_p[i]], 'r*')
            axtrajRZ.plot(self.partdict[ind_p[i]]['R'], self.partdict[ind_p[i]]['z'], 'r-') 
        self._plot_RZsurf(axtrajRZ)
        axtrajRZ.set_xlabel(r'R (m)')
        axtrajRZ.set_ylabel(r'z (m)')
        
        axtrajxy = f.add_subplot(122)        
        for i in ind:
            ll = ind_t[i]
            xi,yi = self.data_i['R'][ll]*np.cos(self.data_i['phi'][ll]), \
                    self.data_i['R'][ll]*np.sin(self.data_i['phi'][ll])
            x,y = self.partdict[ll]['R']*np.cos(self.partdict[ll]['phi']), \
                    self.partdict[ll]['R']*np.sin(self.partdict[ll]['phi'])  
                    
            axtrajxy.plot(xi,yi, 'k*')
            axtrajxy.plot(x,y, 'k-') 
            ll = ind_p[i]
            xi,yi = self.data_i['R'][ll]*np.cos(self.data_i['phi'][ll]), \
                    self.data_i['R'][ll]*np.sin(self.data_i['phi'][ll])
            x,y = self.partdict[ll]['R']*np.cos(self.partdict[ll]['phi']), \
                    self.partdict[ll]['R']*np.sin(self.partdict[ll]['phi'])             
            axtrajxy.plot(xi,yi, 'r*')
            axtrajxy.plot(x,y, 'r-') 

        axtrajxy.set_xlabel(r'x (m)')
        axtrajxy.set_ylabel(r'y (m)')        
        
class h5_particles(particles):
    """
    Class (inherited from particles) handling h5 files (e.g. bbnbi.h5, ascot.h5)
    """
    def __init__(self, fname):
        """
        READ *.h5 FILE
        """
        particles.__init__(self)
        self.fname = fname
        dataf = h5py.File(self.fname)
        self.id = self.fname[-9:-3]
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

        print "FIELDS IN FILE ",self.fname,", section inistate and ",self.endgroup," :"
        print self.field

        # Swapping R and Rprt, since they are opposite to matlab
        self.data_e['R']   = self.data_e['Rprt']
        self.data_e['z']   = self.data_e['zprt']
        self.data_e['phi'] = self.data_e['phiprt']


    def calc_shinethr(self):
        """
        Method to calculate the shine-through with the weights
        """
        
        if self.endgroup != 'shinethr':
            print "WARNING! Check input file, which is ", self.fname

        self._calc_weight()
        self.shinethr=np.zeros((len(self.origins)), dtype=float)
        self.shinethr_abs=np.zeros((len(self.origins)), dtype=float)

        for ind_i,i in enumerate(self.origins):
            ind = self.data_e['origin'][:]==i
            if len(self.data_e['origin'][ind])==0:
                w=0
                e=0
                power = 1
            else:   
                e = self.data_e['energy'][ind][0]
                e = e*1.612e-19
                w = np.sum(self.data_e['weight'][ind_i])
                #wtot = self.w[ind_i]
                power = self.originpower[str(int(i))]*1e6
            self.shinethr[ind_i]=float(e)*w/power
            self.shinethr_abs[ind_i] = float(e)*w
            print "NBI:","{}".format(self.origindict[str(int(i))]),\
                  "Shine-through:", "{0:5f} %".format(self.shinethr[ind_i]*100),\
                  "||  {0:3f} W".format(e*w)


