"""
matteo.vallar@igi.cnr.it - 11/2017

Class for distributions
two classes inside:
distribution_1d(h5 ascot file)
distribution_2d(h5 ascot file with 2D distributions in it)

02/2018: PORTING TO python3 with backwards compatibility using future package
http://python-future.org/index.html

"""
from __future__ import print_function
from utils.fdist_superclass import fdist_superclass

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os.path, math, time
import collections
from utils.plot_utils import _plot_2d, _plot_1d

class distribution_2d(fdist_superclass):

    def __init__(self, infile_n):
        fdist_superclass.__init__(self, infile_n)
        
    def _read(self):
        if os.path.isfile(self.fname) is False:
            print("File ", self.fname, " doesn't exists!")
            raise Exception()

        self.infile=h5py.File(self.fname, 'r')
        if "rzPitchEdist" not in self.infile['distributions'].keys():
            print("No rzPitchE dist in ", self.fname)
        if "rhoPhiPEdist" not in self.infile['distributions'].keys():
            print("No rhoPhiPitchE dist in ", self.fname)
        if "rzMuEdist" not in self.infile['distributions'].keys():
            print("No rzMuE dist in ", self.fname)

        indd = self.fname[-9:-3]
        self.id = indd
        if indd[0:2] != '00':
            self.id = self.fname[-11:-3]
        self._readwall()
        self._RZsurf()

    def plot_space(self, ax=0):
        try:
            self.f_Ep_int.mean()
        except:
            self._integrate_Ep()
       
        z = self.f_Ep_int
        if 'R' in self.dict_dim and 'z' in self.dict_dim:
            x = self.dict_dim['R']
            y = self.dict_dim['z']
            wallrz = [self.R_w, self.z_w]
            surf = [self.Rsurf, self.zsurf, self.RZsurf]
            _plot_2d(x, y, 'R [m]', 'z [m]', Id = self.id, \
                     wallrz=wallrz, surf=surf, ax=ax, dist=z*1e-18, cblabel=r'n ($10^{18}$/$m^3$)')
        elif 'rho' in self.dict_dim and 'phi' in self.dict_dim:
            x = self.dict_dim['rho']
            y = z[0,:]*1.602e-19
            _plot_1d(x,y, r'$\rho$', 'fdist', ax=ax)
            

    def integrate_range_Ep(self, dim_range):
        """
        Integrates over a range of the 2D field
        input: dim_range=[[xmin, xmax], [ymin, ymax]]
        """
        try:
            ftoint = self.f_space_int
        except:
            self._integrate_space()
            ftoint = self.f_space_int
            
        x = self.dict_dim['E']
        y = self.dict_dim['pitch']
        xlim = np.asarray(dim_range)[0,:]
        if min(xlim)>5:
            xlim = xlim*1.6e-19
        ylim = np.asarray(dim_range)[1,:]
        ind_x = np.where((x<=max(xlim)) & (x>=min(xlim)))[0]
        x = x[ind_x][:]
        ind_y = np.where((y<=max(ylim)) & (y>=min(ylim)))[0]
        y = y[ind_y][:]
        int_E = np.trapz(ftoint[ind_x], x, axis = 0)
        int_pitchE = np.trapz(int_E[ind_y], y, axis=0)
        
        print("Fraction of particles in range ",xlim/1.6e-19*1e-6,\
              " MeV :", int_pitchE/self.norm*100., " %")
 

    def _integrate_Ep(self):
        """
        Hidden method to integrate over (E,p)
        """
        dist_toint = self.fdist_notnorm[0,:,:,:,:]/self.norm
#        int_p = np.trapz(dist_toint, self.dict_dim['pitch'], axis=1)        
        int_E = np.trapz(dist_toint, self.dict_dim['E'], axis=0)
        self.f_Ep_int = np.trapz(int_E, self.dict_dim['pitch'], axis=0)
        totvol = np.sum(self.vol)
        self.f_Ep_int=self.f_Ep_int/totvol**2
#        self.f_Ep_int = np.trapz(int_E, self.dict_dim['pitch'], axis=0)
        

        if 'rho' in self.dict_dim.keys():
            for i,el in enumerate(self.vol):
                self.f_Ep_int[0,i] *= el*el 
        # elif 'R' in self.dict_dim and 'z' in self.dict_dim:
        #     dZ=np.abs(self.dict_dim['z'][-1]-self.dict_dim['z'][-2])
        #     dR=np.abs(self.dict_dim['R'][-1]-self.dict_dim['R'][-2])
        #     for i, el in enumerate(self.dict_dim['R']):
        #        self.f_Ep_int[i,:] /= (2*math.pi*el*dZ*dR)  



    def plot_Epitch(self, ax=0, ylim=[0, 600]):
        """
        plot 2D (pitch, energy, int_space(fdist))
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()

        x = self.dict_dim['pitch']
        y = self.dict_dim['E']*1e-3/1.6e-19
        z = self.f_space_int*1.6e-19*1e-17*self.norm
        
        _plot_2d(x, y, dist=z, xlabel=r'$\xi$', ylabel='E [keV]', ax=ax, cblabel=r'$10^{17}$/keV', ylim=ylim)

    def plot_Emu(self):
        """
        plot 2D (mu, energy, int_space(fdist))
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()

        self.xplot = self.dict_dim['mu']/1.6e-19
        self.yplot = self.dict_dim['E']/1.6e-19
        self.zplot = self.f_space_int
        self._plot_2d('mu', 'E', wallrz=0)
        
    def write_pitchE(self):
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        self._write('pitch','E', self.f_space_int, units=['adimensional', 'J'])
            
    def _backup_file(self, fname):
        """
        Does backup of a file fname adding ~
        """
        if os.path.isfile(fname):
            os.rename(fname, fname+'~')
            print("Copy ", fname, " to ", fname+'~ \n')    

    def _write(self, *args, **kwargs):
        """
        Method to write the distribution to dat file
        """
        try:
            units = kwargs['units']
        except:
            units = ['adimensional','adimensional']
        
        x_labels = [args[i] for i in range(len(args)-1)]
        self.y = args[-1]
        fname = self.id+'_'+args[0]+args[1]+'.dat'
        self._backup_file(fname)

        self.info = '' + self.fname + ' ' + ' matteo.vallar@igi.cnr.it ' + \
                        time.strftime("%d/%m/%Y")
        self.info2 = 'For each dimension: name  units   number of bins     min     max'
        self.header = '' 
        for i in range(len(args)-1):
            self.header += args[i]+' '
            self.header += units[i]+' '            
            self.header += ' '+str(len(self.dict_dim[x_labels[i]]))
            self.header += ' '+str(min(self.dict_dim[x_labels[i]]))
            self.header += ' '+str(max(self.dict_dim[x_labels[i]]))
            self.header += ' '
            
        self.header += '\n'
        self.header += "# Normalisation : {0:.5e}".format(round(self.norm, 2))
                    
        with open(fname,'w') as f_handle:
            f_handle.write('# '+self.info+'\n')
            f_handle.write('# '+self.info2+'\n') 
            f_handle.write('# '+self.header+'\n')
            #for lab in x_labels:
            #    np.savetxt(f_handle, self.dict_dim[lab])
            np.savetxt(f_handle, self.y, fmt='%.5e')  

    def _write4d(self, *args, **kwargs):
        """
        Method to write the distribution (4d) to dat file
        """
        try:
            units = kwargs['units']
        except:
            units = ['adimensional','adimensional',
                     'adimensional','adimensional']
        
        x_labels = [args[i] for i in range(len(args)-1)]
        self.y = args[-1]
        fname = self.id+'_'+args[0]+args[1]+'.dat'
        self._backup_file(fname)

        self.info = '' + self.fname + ' ' + ' matteo.vallar@igi.cnr.it ' + \
                        time.strftime("%d/%m/%Y")
        self.info2 = 'For each dimension: name  units   number of bins     min     max'
        self.header = '' 
        for i in range(len(args)-1):
            self.header += args[i]+' '
            self.header += units[i]+' '            
            self.header += ' '+str(len(self.dict_dim[x_labels[i]]))
            self.header += ' '+str(min(self.dict_dim[x_labels[i]]))
            self.header += ' '+str(max(self.dict_dim[x_labels[i]]))
            self.header += ' '
            
        self.header += '\n'
        self.header += "# Normalisation : {0:.5e}".format(round(self.norm, 2))
                    
        f_handle = open(fname,'w')
        f_handle.write('# '+self.info+'\n')
        f_handle.write('# '+self.info2+'\n') 
        f_handle.write('# '+self.header+'\n')
        #for lab in x_labels:
        #    np.savetxt(f_handle, self.dict_dim[lab])
        for sl in self.y[0,:,:,:,:]:
            for ssll in sl:
                np.savetxt(f_handle, ssll, fmt='%.5e')
                f_handle.write('\n')
            f_handle.write('\n')

    def _computenorm(self):
        """
        calculates the norm of the function
        """
        if "pitch" in self.dict_dim:
            try:
                self.f_spacep_int()
            except:
                self._integrate_spacep()
                
            self.norm = np.trapz(self.f_spacep_int, self.dict_dim['E'])
            print("NORM = ", self.norm)
            #self.norm=1
        elif "mu" in self.dict_dim:
            try:
                self.f_spacemu_int()
            except:
                self._integrate_spacemu()
            #print "NORM = ", self.norm        

    def build_fdist(self):
        """
        Method to read the ordinates of the 2d distribution
        """
        try:
            self.dict_dim
        except:
            print("No dictionary of dimensions created")
            raise ValueError

        self.norm = 1
        # 6th dimension is the one labelling the beams
        tmp = self.dist_h5['ordinate'].value
#        print(np.shape(tmp))
#        for i in range(np.shape(tmp)[0]):
#            f=plt.figure()
#            ax=f.add_subplot(111)
#            distplot=np.sum(tmp, axis=-2)
#            distplot=np.sum(distplot, axis=-2)
#            f.suptitle(str(i))
#            ax.contourf(distplot[i,0,:,:,0])
#        tmp = tmp[[6,9,11,10],:,:,:,:,:,:]    
        fdist = np.sum(tmp, axis=0)[:,:,:,:,:,0] #time, E,pitch,z,r
        self.fdist_notnorm = fdist
        self._computenorm()
        #self.fdist_norm = self.fdist_notnorm/self.norm

    def collect_dim(self):
        """
        methods to read the dictionary of the dimensions of the 4D distributions and store it
        """
        self._read_dim_h5()
        self._fix_dim()

    def _fix_dim(self):
        """
        Hidden method to make the abscissae the correct length (same as ordinate)
        """
        try:
            self.dict_dim
        except:
            self._read_dim()

        for dim in self.dict_dim:
            tmp_dim = self.dict_dim[dim]
            self.dict_dim[dim] = np.linspace(min(tmp_dim), max(tmp_dim), len(tmp_dim)-1)

    def _RZsurf(self):
        """
        Reads the position of RZ surfaces from ascot file
        now the edge is set to the value for scenario 5 from JT60SA
        """
        f = self.infile
        self.RZsurf = f['bfield/2d/psi'].value
        self.Rsurf = f['bfield/r']
        self.zsurf = f['bfield/z']
        edge = f['boozer/psiSepa'][:]; axis=f['boozer/psiAxis'][:]
        self.RZsurf = (-1*self.RZsurf - axis )/(edge-axis)
        self.RZsurf = np.sqrt(self.RZsurf)            
        
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
        
class frzpe(distribution_2d):
    
    def __init__(self, infile_n):
        """
        Module to initialise the distributions (now works only with rzPitchEdist)
        self.fdist['ordinate'] has the following shape: (ind beam, time, energy, pitch, z, R, #ions)
        """
        distribution_2d.__init__(self, infile_n)
        try:
            self.dist_h5 = self.infile['distributions/rzPitchEdist'] 
        except:
            raise ValueError
        self.__name__ = 'frzpe'

        self.dict_dim = collections.OrderedDict([('R',[]),('z',[]),('pitch',[]),('E',[]),('t',[])])
        self.collect_dim()
        self.build_fdist()
        self._integrate_space()        

    def _integrate_space(self):
        """
        Function to integrate over (R,z)
        """
        dist_toint = self.fdist_notnorm[0,:,:,:,:]/self.norm

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #E,pitch

    
    def plot_space_enslice(self, sliceind):
        """
        Function to plot over (R,z) on a E defined
        """
        self._integrate_pitch_enslice(sliceind)
        self.zplot = self.f_xi_int
        if 'R' in self.dict_dim and 'z' in self.dict_dim:
            self.xplot = self.dict_dim['R']
            self.yplot = self.dict_dim['z'] 
            _plot_2d('R [m]', 'z [m]', wallrz=1, surf=1, \
                          title=str(sliceind))

    def _integrate_pitch_enslice(self, sliceind):
        """
        Function to integrate over (R,z) on a E defined
        """
        sliceind*=1.602e-19
        ind_E = np.argmin(self.dict_dim['E']-sliceind < 0)
        print(ind_E, self.dict_dim['E'][ind_E])
        dist_toint = self.fdist_notnorm[0,ind_E,:,:,:]/self.norm

        int_xi  = np.trapz(dist_toint, self.dict_dim['pitch'], axis = 0)
        self.f_xi_int = int_xi #E,pitch
        
        """
        Function to integrate over (R,z) on a E defined
        """
        sliceind/=1.602e-19
        ind_E = np.argmin(self.dict_dim['E']-sliceind < 0)
        dist_toint = self.fdist_notnorm[0,ind_E,:,:,:]/self.norm

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #E,pitch
        
    def _integrate_space_pslice(self, sliceind):
        """
        Function to integrate over (R,z) on a pitch defined
        """
        dist_toint = self.fdist_notnorm[0,:,sliceind,:,:]/self.norm

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #E,pitch

    def plot_RZposition(self, sliceR, slicez, **kwargs):
        """
        makes a plot of E, pitch on a R,z position, given in sliceR and sliceZ
        """
        ind_R = np.argmin(self.dict_dim['R']-sliceR < 0)
        print(ind_R, self.dict_dim['R'][ind_R])
        ind_z = np.argmin(self.dict_dim['z']-slicez < 0)
        dist_toplot = self.fdist_notnorm[0,:,:,ind_R, ind_z]/self.norm

        self.xplot = self.dict_dim['pitch']
        self.yplot = self.dict_dim['E']/1.6e-19*1e-3
        self.zplot = dist_toplot
        if 'fname' in kwargs:
            self._plot_2d(r'$\xi$', 'E [keV]', \
                          title='R='+str(sliceR)+' z='+str(slicez), \
                          fname=kwargs['fname'])
        else:
            self._plot_2d(r'$\xi$', 'E [keV]', \
                          title='R='+str(sliceR)+' z='+str(slicez))

    def _get_Eposition(self, sliceE):
        """
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        ind_E = np.argmin(self.dict_dim['E']/1.602e-19-sliceE*1000. < 0)
        dist_toplot = self.f_space_int[ind_E, :]        
        return dist_toplot
    
    def plot_Eposition(self, sliceE, **kwargs):
        """
        makes a plot of a slice of energy (sliceE is in keV) vs pitch integrated over 
        configuration space
        """
        dist_toplot = self._get_Eposition(sliceE)

        self.xplot = self.dict_dim['pitch']
        self.yplot = dist_toplot
        if 'ax' in kwargs:
            self._plot_1d(r'$\xi$', r'E [keV]',\
                          ax=kwargs['ax'])            
        if 'fname' in kwargs:
            self._plot_1d(r'$\xi$', r'E [keV]',\
                          title='E='+str(sliceE), \
                          fname=kwargs['fname'])
        else:
            self._plot_1d(r'$\xi$',r'E [keV]', title='E='+str(sliceE))        

    def write_allf(self):
        self._write4d('R','z','pitch','E', self.fdist_notnorm/self.norm, \
                         units=['m','m','adimensional', 'J'])
    def write_allf_notnorm(self):
        self._write4d('R','z','pitch','E', self.fdist_notnorm, \
                       units=['m','m','adimensional', 'J'])

            
class frhophipe(distribution_2d):
    
    def __init__(self, infile_n):
        """
        Module to initialise the distributions (now works only with rzPitchEdist, rhophipitchE)
        self.fdist['ordinate'] has the following shape: (ind beam, time, energy, pitch, phi, rho, #ions)
        """
        distribution_2d.__init__(self, infile_n)
        try:
            self.dist_h5 = self.infile['distributions/rhoPhiPEdist'] 
        except:
            raise ValueError
        self.__name__ = 'frhophipe'
        
        self.dict_dim = collections.OrderedDict([('rho',[]),('phi',[]),('pitch',[]),('E',[]),('t',[])])
        self.vol = self.infile['distributions/rhoDist/shellVolume'].value

        self.collect_dim()
        self.build_fdist()
        self._integrate_space()

    def _integrate_space(self):
        """
        Function to integrate over (rho,phi)
        """
        dist_toint = self.fdist_notnorm[0,:,:,:,:]/self.norm

        #np.cumsum(shellVol) is the profile of the volume, enclosed in a rho surf
        for i,el in enumerate(np.cumsum(self.vol)):
            dist_toint[:,:,:,i] *= el/self.shape_dim['phi']       
            
        int_rho    = np.trapz(dist_toint, self.dict_dim['rho'], axis = -1)
        
        #int_rhophi  = np.trapz(int_rho   , self.dict_dim['phi'], axis = -1) 
        int_rhophi = int_rho[:,:,0]*2.*np.pi
        self.f_space_int = int_rhophi #E,pitch  
        
    def write_allf(self):
        self._write4d('rho','phi','pitch','E', self.fdist_notnorm/self.norm, \
                         units=['adim','rad','adim', 'J'])

    def write_allf_notnorm(self):
        self._write4d('rho','phi','pitch','E', self.fdist_notnorm, \
                         units=['adim','rad','adim', 'J'])

    def plot_rhoposition(self, slicerho,**kwargs):
        """
        makes a plot of E, pitch on a rho position, given in slicerho
        """
        ind_rho = np.argmin(self.dict_dim['rho']-slicerho< 0)
        print(ind_rho, self.dict_dim['rho'][ind_rho])
        dist_toplot = self.fdist_notnorm[0,:,:,0,ind_rho]/self.norm

        self.xplot = self.dict_dim['pitch']
        self.yplot = self.dict_dim['E']/1.6e-19*1e-3
        self.zplot = dist_toplot
        if 'fname' in kwargs:
            _plot_2d(r'$\xi$', 'E [keV]', \
                          title=r'$\rho$='+str(slicerho), \
                          fname=kwargs['fname'])
        else:
            _plot_2d(r'$\xi$', 'E [keV]', title=r'$\rho$='+str(slicerho))

       
class frzv(distribution_2d):
    
    def __init__(self, infile_n):
        """
        Module to initialise the distributions RZMUEDIST
        self.fdist['ordinate'] has the following shape: (ind beam, time, energy, pitch, z, R, #ions)
        """
        distribution_2d.__init__(self, infile_n)
        try:
            self.dist_h5 = self.infile['distributions/rzVDist'] 
        except:
            raise ValueError
        self.__name__ = 'frzv'

        self.dict_dim = collections.OrderedDict([('R',[]),('z',[]),('vpara',[]),('vperp',[]),('t',[])])
        self.collect_dim()
        self.build_fdist()
        self._integrate_space()        
        
    def _integrate_space(self):
        """
        Function to integrate over (R,z)
        """
        dist_toint = self.fdist_notnorm[0,:,:,:,:]/self.norm

        for i, el in enumerate(self.dict_dim['R']):
            dist_toint[:,:,:,i] *= 2*math.pi*el

        int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
        int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
        self.f_space_int = int_Rz #E,pitch

    def _integrate_v(self):
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()

        dist_toint = self.f_space_int
        int_vpara = np.trapz(dist_toint, self.dict_dim['vpara'], axis=0)
        self.norm = np.trapz(int_vpara, self.dict_dim['vperp'], axis=0)
        
        
        
    def plot_vv(self, ax=0):
        """
        plot 2D (vpara, vperp, int_space(fdist))
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()

        x = self.dict_dim['vpara']
        y = self.dict_dim['vperp']
        z = self.f_space_int*self.norm
        
        _plot_2d(x*1e-6, y*1e-6, dist=z, xlabel=r'v$_{\parallel}$[km/ms]', ylabel=r'v$_{\perp}$[km/ms]', ax=ax)
        plt.axis('equal')
