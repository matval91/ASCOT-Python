"""
matteo.vallar@igi.cnr.it
Classes to handle magnetic fields I/O for ascot&bbnbi

"""
from __future__ import print_function

import numpy as np
import h5py, math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import ReadEQDSK
from scipy.interpolate import griddata
import scipy.optimize
import scipy.interpolate as interp
import scipy.io as sio
from ascot_utils import common_style

class Bfield_ascot:
    """ Handles Bfield store in h5 file

    Class for handling the magnetic field specifications and plots
    Better to ignore boozer field, it is useless

    Parameters:
        infile_name (str): name of file to analize
    Attributes:
        None
    Notes:    
    |  Groups in h5 file (09/01/2017, ascot4)
    |  /bfield                  Group
    |  /bfield/2d               Group
    |  /bfield/2d/bphi          Dataset {600, 300}
    |  /bfield/2d/psi           Dataset {600, 300}
    |  /bfield/nphi             Dataset {1}
    |  /bfield/nr               Dataset {1}
    |  /bfield/nz               Dataset {1}
    |  /bfield/r                Dataset {300}
    |  /bfield/raxis            Dataset {1}
    |  /bfield/z                Dataset {600}
    |  /bfield/zaxis            Dataset {1}


    |  /boozer                  Group
    |  /boozer/Babs             Dataset {200, 100}
    |  /boozer/Ifunc            Dataset {100}
    |  /boozer/R                Dataset {200, 100}
    |  /boozer/axisR            Dataset {1}
    |  /boozer/axisz            Dataset {1}
    |  /boozer/btheta           Dataset {200, 100}
    |  /boozer/delta            Dataset {200, 100}
    |  /boozer/gfunc            Dataset {100}
    |  /boozer/nu               Dataset {200, 100}
    |  /boozer/psi              Dataset {100}
    |  /boozer/psiAxis          Dataset {1}
    |  /boozer/psiMax           Dataset {1}
    |  /boozer/psiMin           Dataset {1}
    |  /boozer/psiSepa          Dataset {1}
    |  /boozer/qprof            Dataset {100}
    |  /boozer/rgridmax         Dataset {1}
    |  /boozer/rgridmin         Dataset {1}
    |  /boozer/theta            Dataset {200}
    |  /boozer/z                Dataset {200, 100}
    |  /boozer/zgridmax         Dataset {1}
    |  /boozer/zgridmin         Dataset {1}

    
    METHODS:
    __init__(self, infile) to store variables from h5 file
    checkplot(self) to plot magnetic values and check them

    HIDDEN METHODS:
    _read_wall_h5(self): stores wall data from h5 file
    _sanitycheck(self): checks the input equilibrium
    """

    
    def __init__(self, infile_name):
        #defining dictionary for handling data from h5file
        self.labdict = {"bphi":"/bfield/2d/bphi","psi_2D":"/bfield/2d/psi",\
                        "nr":"/bfield/nr", "nz":"/bfield/nz",\
                        "r":"/bfield/r", "z":"/bfield/z",\
                        "raxis":"/bfield/raxis", "zaxis":"/bfield/zaxis",\
                        "q":"/boozer/qprof", "psi_rho":"/boozer/psi"}
        
#                         "Babs":"/boozer/Babs", "I":"/boozer/Ifunc",\
#                         "r_boo":"/boozer/R","axisR_boo":"/boozer/axisR", "axisZ_boo":"/boozer/axisz",\
#                         "btheta":"/boozer/btheta", "delta":"/boozer/delta",\
#                         "g":"/boozer/gfunc", "nu":"/boozer/nu",\
#                         "psi_rho":"/boozer/psi", "psi_ax":"/boozer/psiAxis",\
#                         "maxPsi":"/boozer/psiMax", "minPsi":"/boozer/psiMin",\
#                         "psi_sepa":"/boozer/psiSepa",\
#                         "rmax":"/boozer/rgridmax", "rmin":"/boozer/rgridmin",\
#                         "theta":"/boozer/theta",\
#                         "z_boo":"/boozer/z","zmax":"/boozer/zgridmax", "zmin":"/boozer/zgridmin"\
#                         }

        self.infile=h5py.File(infile_name)
        self.infile_n = infile_name
        self.vardict = {}
        #store data with the correct labels
        
        for k in self.labdict:
            self.vardict[k] = self.infile[self.labdict[k]]

        self.nrho = self.vardict['psi_rho'].shape[0]
        self.rho = self.vardict['psi_rho'][:]**0.5

        #these nR and nZ are from the group /bfield/, the ones for /boozer/ are on rho and theta (spacing of grid)
        self.nR = self.vardict['nr'][:]        
        self.nZ = self.vardict['nz'][:]
        self.R = self.vardict['r'][:]
        self.Z = self.vardict['z'][:]
        self.rmin=np.min(self.R)
        self.rmax=np.max(self.R)
        self.zmin=np.min(self.Z)
        self.zmax=np.max(self.Z)

        self._read_wall_h5()

    def _read_wall_h5(self):
        """stores wall data from h5 file
        
        Reads wall data from ascot.h5 file
        
        Parameters:
            None
        Attributes:
            None
        Note:
            Could be implemented in ascot_utils

        """
        self.walllabdict = {"R_wall":"/wall/2d/R", "Z_wall":"/wall/2d/z",\
                            "divflag":"/wall/2d/divFlag", "segLen":"/wall/2d/segLen"}
        self.w = dict.fromkeys(self.walllabdict)
        for k in self.walllabdict:
            self.w[k] = self.infile[self.walllabdict[k]].value
        self.R_w = self.w['R_wall']
        self.z_w = self.w['Z_wall']

    def _sanitycheck(self):
        """ checks the input equilibrium
        
        In the group sanity check the magnetic field is stored (just the values needed
        for ASCOT to work), so here you can take the data to check you're doing things
        right
        
        Parameters:
            None
        Attributes:
            None

        """
        self.sanitydict = {'Ip': None, 'b0': None,\
                           'bphi': None,'br': None,'bz': None,'psiAtAxis': None,\
                           'psiAtSeparatrix': None,'rmax': None, 'rmin': None,\
                           'separatrix': None,'zmax': None,'zmin': None}
        for key in self.sanitydict:
            self.sanitydict[key] = self.infile['sanityChecks/'+key].value
        for key in ['Ip', 'b0','psiAtAxis', 'psiAtSeparatrix',\
                    'rmin','rmax','zmin','zmax']:
            print(key, " ", self.sanitydict[key])
            

    def _plot_2dpsi(self,ax):
        """
        
        """
        edge = self.infile['boozer/psiSepa'][:]; axis=self.infile['boozer/psiAxis'][:]
        edge = self.sanitydict['psiAtSeparatrix']; axis=self.sanitydict['psiAtAxis']

        yplot = self.vardict['psi_2D'][:]
        r,z = self.vardict['r'], self.vardict['z']

        #ax.set_title('Poloidal flux')
        CS = ax.contour(r,z,yplot, 30)
        ax.contour(r,z,yplot,[-edge], colors='k', linewidths=2.5, linestyles='--')
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        #CB = plt.colorbar(CS)
        if self.R_w[0]!=0:
            ax.plot(self.R_w, self.z_w, 'k',linewidth=2.5)
        ax.axis('equal')
        ax.grid('on')
        plt.setp(ax.get_yticklabels()[0], visible=False) 
        M = 4
        yticks = ticker.MaxNLocator(M)
        xticks = ticker.MaxNLocator(M)
        # tick positions with set_minor_locator.
        ax.yaxis.set_major_locator(yticks)
        #ax.yaxis.set_minor_locator(yticks_m)
        ax.xaxis.set_major_locator(xticks)

    def checkplot(self, f=0):
        """Plots psi_2D, q, Bphi
        
        Method to plot the values and check the magnetic field we are looking at

        Parameters:
            f (obj): figure object where to plot
        Attributes:
            None
        
        """
        try:
            self.sanitydict['psiAtSeparatrix']
        except:
            self._sanitycheck()

        common_style()

        if f==0:
            f = plt.figure(figsize=(20, 8))
            #f.text(0.01, 0.01, self.infile_n)
            ax2d = f.add_subplot(131)
            axq = f.add_subplot(132)
            axf = f.add_subplot(133)
        else:
            f=plt.gcf()
            ax2d = f.axes[0]
            axq  = f.axes[2]
            axf  = f.axes[3]
            plt.rc('lines', linestyle='--')
        ax2d = f.add_subplot(131)

        self._plot_2dpsi(ax2d)

        M = 4
        yticks = ticker.MaxNLocator(M)
        xticks = ticker.MaxNLocator(M)
        # tick positions with set_minor_locator.
        ax2d.yaxis.set_major_locator(yticks)
        #ax.yaxis.set_minor_locator(yticks_m)
        ax2d.xaxis.set_major_locator(xticks)

        axq = f.add_subplot(132)
        axq.plot(self.rho, np.abs(self.vardict['q']), lw=2.3, color='k')
        axq.set_xlabel(r'$\rho_{POL}$')
        axq.set_ylabel(r'q'); axq.set_ylim([0,6])
        axq.grid('on')
        M = 4
        yticks = ticker.MaxNLocator(M)
        xticks = ticker.MaxNLocator(M)
        # tick positions with set_minor_locator.
        axq.yaxis.set_major_locator(yticks)
        #ax.yaxis.set_minor_locator(yticks_m)
        axq.xaxis.set_major_locator(xticks)
        axf = f.add_subplot(133)
        axf.plot(self.vardict['r'],  self.vardict['bphi'][50,:], lw=2.3, color='k')
        axf.set_xlabel(r'R [m]')
        axf.set_ylabel(r'$B_\phi$')
        axf.grid('on')
        f.tight_layout()
        plt.rc('lines', linestyle='-')

        plt.setp(axq.get_yticklabels()[0], visible=False) 
        plt.setp(axf.get_yticklabels()[0], visible=False) 
        M = 4
        yticks = ticker.MaxNLocator(M)
        xticks = ticker.MaxNLocator(M)
        # tick positions with set_minor_locator.
        axf.yaxis.set_major_locator(yticks)
        #ax.yaxis.set_minor_locator(yticks_m)
        axf.xaxis.set_major_locator(xticks)

        plt.show()

        
class Bfield_eqdsk:
    """ Class handling eqdsk magnetic field

    Script for writing and reading the magnetic background
    porting from matlab (16-1-2017)

    Parameters:
        |  infile (str): filename of the eqdsk (with only 4 stripped strings in the first line before the nR and nZ)
        |  nR (int):  number of R grid to output. Usually 259
        |  nz (int):  number of z grid to output. Usually 259
        |  devnam (str): name of the device (JT60SA, TCV)
        |  COCOS (int): number identifying COCOS.
    
    Attributes:
        None
    
    Methods:
        |  eqdsk_checkplot : Method to plot the values (2D psi, q, poloidal flux) and check the magnetic field we are looking at
        |  write: Function calling the two methods to write the header and the bkg
        |  write_bkg: Write to input.magn_bkg file
        |  write_head: Write to input.magn_header file
        |
        |  build_lim: Function calling the two methods to build the header (with limiter) and bkg dictionary
        |  build_SN: Function calling the two methods to build the header (with SN) and bkg dictionary
        |  build_header_lim: Method to build header file from eqdsk without single nulls (one point in PFxx) 
        |  build_header_SN: Method to build header file from eqdsk with one single null (two points in PFxx)
        |  
        |  calc_field: Function to calculate toroidal fields (fields on poloidal plane set to 0
    """
    
    def __init__(self, infile, nR, nz, devnam, COCOS, *args):
        self.COCOS = COCOS
        self.devnam = devnam
        self._read_wall()
        self._import_from_eqdsk(infile)
        #these are the dimensions of the output arrays
        self.nR=nR
        self.nz=nz
        if len(args)!=0:
            self.dr, self.dz = np.asarray(args[0]), np.asarray(args[1])
            print("Shifting of quantity ", self.dr, " and ", self.dz)
            self._shift_eq()


    def _import_from_eqdsk(self, infile_eqdsk):
        """ importing from eqdsk
        function for import from eqdsk file

        Parameters:
            infile_eqdsk (str): name of eqdsk file to read
        Attributes:
            None
        Notes:
        these are the data of the eqdsk struct:
        
            self.comment=comment
            self.switch=switch
            self.nrbox=nrbox
            self.nzbox=nzbox
            self.rboxlength=rboxlength
            self.zboxlength=zboxlength
            self.R0EXP=R0EXP
            self.rboxleft=rboxleft
            self.Raxis=Raxis
            self.Zaxis=Zaxis
            self.psiaxis=psiaxis
            self.psiedge=psiedge
            self.B0EXP=B0EXP
            self.Ip=Ip
            self.T=T
            self.p=p
            self.TTprime=TTprime
            self.pprime=pprime
            self.psi=psi
            self.q=q
            self.nLCFS=nLCFS
            self.nlimits=nlimits
            self.R=R
            self.Z=Z
            self.R_limits=R_limits
            self.Z_limits=Z_limits
            self.R_grid=R_grid
            self.Z_grid=Z_grid
            self.psi_grid=psi_grid
            self.rhopsi=rhopsi
        
        """    
        self.eqdsk= ReadEQDSK.ReadEQDSK(infile_eqdsk)
        self.infile = infile_eqdsk
        self.eqdsk.psi = np.reshape(self.eqdsk.psi, (self.eqdsk.nzbox, self.eqdsk.nrbox))       
        self.R_eqd = np.linspace(self.eqdsk.rboxleft, self.eqdsk.rboxleft+self.eqdsk.rboxlength, self.eqdsk.nrbox)
        self.Z_eqd = np.linspace(-self.eqdsk.zboxlength/2., self.eqdsk.zboxlength/2., self.eqdsk.nzbox)
        self.cocos_transform(self.COCOS)
        #This is for Equilibrium from CREATE for scenario 2, also to modify in build bkg
        self.psi_coeff = interp.RectBivariateSpline(self.Z_eqd, self.R_eqd, self.eqdsk.psi)
        #self.psi_coeff = interp.RectBivariateSpline(self.R_eqd, self.Z_eqd, self.eqdsk.psi) 
       
    def cocos_transform(self, COCOS):
        """ cocos transformations
        This function converts the magnetic input from their starting cocos to eqdsk 5 (needed by ascot).

        Parameters:
            COCOS (int): input cocos. Now useable only 2,3,4,7,12,13,14,17
        Attributes:
            None
        
        """

        print("COCOS tranformation from "+str(COCOS)+" to 5")
        cocos_keys = ['sigma_Bp', 'sigma_RphiZ', 'sigma_rhothetaphi', 'sign_q_pos', 'sign_pprime_pos', 'exp_Bp']
        pi = math.pi
        cocosin = dict.fromkeys(cocos_keys)
        cocosin['exp_Bp'] = 0 # if COCOS>=10, this should be 1
        if COCOS>10:
            cocosin['exp_Bp'] = +1 # if COCOS>=10, this should be 1
            
        
        if COCOS==2 or COCOS==12:
            #These cocos are for CHEASE (2) - IN
            cocosin['sigma_Bp'] = +1
            cocosin['sigma_RphiZ'] = -1
            cocosin['sigma_rhothetaphi'] = +1
            cocosin['sign_q_pos'] = +1
            cocosin['sign_pprime_pos'] = -1
        elif COCOS==3 or COCOS==13:
            #These cocos are for EFIT (3) - IN
            cocosin['sigma_Bp'] = -1
            cocosin['sigma_RphiZ'] = +1
            cocosin['sigma_rhothetaphi'] = -1
            cocosin['sign_q_pos'] = -1
            cocosin['sign_pprime_pos'] = +1
        elif COCOS==4 or COCOS==14:
            cocosin['sigma_Bp'] = -1
            cocosin['sigma_RphiZ'] = -1
            cocosin['sigma_rhothetaphi'] = -1
            cocosin['sign_q_pos'] = -1
            cocosin['sign_pprime_pos'] = +1
        elif COCOS==7 or COCOS==17:
            #These cocos are for LIUQE(17) - IN
            cocosin['sigma_Bp'] = -1
            cocosin['sigma_RphiZ'] = +1
            cocosin['sigma_rhothetaphi'] = +1
            cocosin['sign_q_pos'] = +1
            cocosin['sign_pprime_pos'] = +1
       
        else:
            print(str(COCOS)+" Not Implemented \n")
            raise ValueError
        
        cocosin['sigma_ip'] = np.sign(self.eqdsk.Ip)
        cocosin['sigma_b0'] = np.sign(self.eqdsk.B0EXP)

        #These cocos are for ASCOT - OUT
        cocosout = dict.fromkeys(cocos_keys)
        cocosout['sigma_Bp'] = 1
        cocosout['sigma_RphiZ'] = +1
        cocosout['sigma_rhothetaphi'] = +1
        cocosout['sign_q_pos'] = +1
        cocosout['sign_pprime_pos'] = -1
        cocosout['exp_Bp'] = 0
        cocosout['sigma_ip'] = +1
        cocosout['sigma_b0'] = -1
        
        # Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
        #sigma_Ip_eff = cocosin['sigma_RphiZ'] * cocosout['sigma_RphiZ']
        #sigma_B0_eff = cocosin['sigma_RphiZ'] * cocosout['sigma_RphiZ']
        # Since we want sigmaip and sigmab0 defined, we must use
        sigma_Ip_eff = cocosin['sigma_ip']*cocosout['sigma_ip']
        sigma_B0_eff = cocosin['sigma_b0']*cocosout['sigma_b0']
        sigma_Bp_eff = cocosin['sigma_Bp'] * cocosout['sigma_Bp']
        exp_Bp_eff = cocosout['exp_Bp'] - cocosin['exp_Bp']
        sigma_rhothetaphi_eff = cocosin['sigma_rhothetaphi'] * cocosout['sigma_rhothetaphi']
        # Define input
        F_in = self.eqdsk.T
        FFprime_in = self.eqdsk.TTprime
        pprime_in = self.eqdsk.pprime
        psirz_in = self.eqdsk.psi
        psiaxis_in = self.eqdsk.psiaxis
        psiedge_in = self.eqdsk.psiedge
        q_in = self.eqdsk.q
        b0_in = self.eqdsk.B0EXP
        ip_in = self.eqdsk.Ip
        
        # Transform
        F = F_in * sigma_B0_eff
        FFprime = FFprime_in*sigma_Ip_eff*sigma_Bp_eff/(2*pi)**exp_Bp_eff
        pprime = pprime_in * sigma_Ip_eff * sigma_Bp_eff / (2*pi)**exp_Bp_eff
        _fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2*pi)**exp_Bp_eff
        psirz = psirz_in * _fact_psi       
        psiaxis = psiaxis_in * _fact_psi
        psiedge = psiedge_in * _fact_psi
        q = q_in * sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff
        b0 = b0_in * sigma_B0_eff
        ip = ip_in * sigma_Ip_eff
        # Define output
        self.eqdsk.T = F
        self.eqdsk.TTprime = FFprime
        self.eqdsk.pprime = pprime
        self.eqdsk.psi = psirz
        self.eqdsk.psiaxis = psiaxis
        self.eqdsk.psiedge = psiedge
        self.eqdsk.q = q
        self.eqdsk.B0EXP = b0
        self.eqdsk.Ip = ip

    def _shift_eq(self):
        """
        NOT FULLY FUNCTIONAL
        Shift the equilibrium of the quantity dr and dz (Dr>0: shift to right, Dz>0: shift up)
        """
        self.eqdsk.R0EXP    += self.dr
        self.eqdsk.rboxleft += self.dr
        self.eqdsk.Raxis    += self.dr
        self.eqdsk.R        += self.dr
        self.eqdsk.R_grid   += self.dr
        self.R_eqd          += self.dr

        self.eqdsk.Zaxis    += self.dz
        self.eqdsk.Z        += self.dz
        self.eqdsk.Z_grid   += self.dz
        self.Z_eqd          += self.dz

    def eqdsk_checkplot(self, f=0):
        """plot of 2D psi, q, bphi

        Method to plot the values (2D psi, q, bphi) and check the magnetic field we are looking at
        
        Parameters:
            f (obj): figure object where to plot. if undefined, f=0
        Attributes:
            None

        """
        try:
            self.param_bphi
        except:
            self.calc_field()
            
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        plt.rc('axes', labelsize=15)
        plt.rc('figure', facecolor='white')
        if f==0:
            f = plt.figure(figsize=(20, 8))
            f.text(0.01, 0.01, self.infile)
        else:
            f=plt.gcf()
        ax2d = f.add_subplot(131)
        r,z = self.R_eqd, self.Z_eqd
        #r = np.linspace(float(np.around(np.min(self.R_w), decimals=2)), float(np.around(np.max(self.R_w), decimals=2)), self.nR)
        #z = np.linspace(float(np.around(np.min(self.z_w), decimals=2)), float(np.around(np.max(self.z_w), decimals=2)), self.nz)

 
        CS = ax2d.contour(r,z, self.psi_coeff(z,r), 30)
        plt.contour(r,z, self.psi_coeff(z,r), [self.eqdsk.psiedge], colors='k', linewidths=3.)

        ax2d.set_xlabel("R [m]")
        ax2d.set_ylabel("Z [m]")
        CB = plt.colorbar(CS)
        if self.R_w[0]!=0:
            ax2d.plot(self.R_w, self.z_w, 'k',linewidth=2)
        ax2d.axis('equal')
        axq = f.add_subplot(132)
        axq.plot(self.eqdsk.rhopsi, self.eqdsk.q, lw=2.3, color='k')
        axq.set_xlabel(r'$\rho_{POL}$')
        axq.set_ylabel(r'q')

        axf = f.add_subplot(133)
        #axf.plot(self.R_eqd, self.eqdsk.T)
        axf.plot(r, self.param_bphi(r,z)[len(r)/2,:], lw=2.3, color='k')
        axf.set_xlabel(r'R [m]')
        axf.set_ylabel(r'$B_\phi$')
        f.tight_layout()
        plt.show()

    def plot_Bfield(self):
        """
        Plots the magnetic field used (only thing ASCOT is interested by)
        """
        try:
            self.Br_t.mean()
        except:
            self.calc_field()

        self.Br_t, self.Bz_t, self.Bphi

        f=plt.figure()
        axr=f.add_subplot(131)
        CS = plt.contour(self.R_eqd, self.Z_eqd, self.Br_t)
        CB = plt.colorbar(CS)
        axr.set_title("R")

        axz = f.add_subplot(132)
        CS = plt.contour(self.R_eqd, self.Z_eqd, self.Bz_t)
        CB = plt.colorbar(CS)
        axz.set_title("z")

        axphi = f.add_subplot(133)
        CS = plt.contour(self.R_eqd, self.Z_eqd, self.param_bphi(self.R_eqd, self.Z_eqd))
        CB = plt.colorbar(CS)
        axphi.set_title("phi")

        plt.show()

    def write(self):
        """
        Function calling the two methods to write the header and the bkg
        """
        self.write_head()
        self.write_bkg()

    def build_lim(self):
        """ limiter building

        Function calling the two methods to build the header (with limiter) and bkg dictionary
        
        Parameters: 
            None
        Attributes:
            None

        """
        try: 
            self.hdr['Vol'].mean()
        except:
            self.build_header_lim()
        
        try:
            self.bkg['Bphi'].mean()
        except:
            self.build_bkg()
  

    def build_SN(self):
        """ building single null

        Function calling the two methods to build the header (with SN) and bkg dictionary
        In this case there are two special points (and the x-point can be found with ginput from plot)

        Parameters:
            None
        Attributes:
            None

        """
        try: 
            self.hdr['Vol'].mean()
        except:
            self.build_header_SN()
        
        try:
            self.bkg['Bphi'].mean()
        except:
            self.build_bkg()
        
    def build_header_lim(self):
        """ building limiter header

        |  Method to build header file from eqdsk without single nulls (one point in PFxx) 
        |  -The first five values of the eqdsk (nShot, tShot, modflg, FPPkat, IpiFPP) are already set correctly  
        |  -The last quantities (rhoPF, PFL, Vol, Area, Qpl) are already set

        Parameters:
            None
        Attributes:
            None
        
        """
        print("Build hdr (limiter)")

        
        nrho = len(self.eqdsk.rhopsi)
        dummy=np.linspace(0,1,nrho)
        
        self.hdr={'nSHOT':0,'tSHOT':0,'modflg':0,'FPPkat':0,'IpiFPP':self.eqdsk.Ip,\
                  'PFxx':[],'RPFx':[],'zPFx':[],'SSQ':[], 'devnam':self.devnam,\
                  'rhoPF':nrho,'PFL':dummy,'Vol':dummy,'Area':dummy,'Qpl':dummy} 
        
        # find axis
        self.ax = self._min_grad(x0=[self.eqdsk.Raxis, self.eqdsk.Zaxis])     
        self.axflux = self.psi_coeff(self.ax[0], self.ax[1])*(2*math.pi)
        print("remember: I am multiplying psi axis times 2pi since in ascot it divides by it!")

        # poloidal flux of the special points (only one in this case)
        self.hdr['PFxx'] = [self.axflux[0][0]]
        self.hdr['RPFx'] = [self.ax[0]]
        self.hdr['zPFx'] = [self.ax[1]]
        self.hdr['SSQ']  = [self.eqdsk.Raxis, self.eqdsk.Zaxis, 0, 0]

        
    def build_header_SN(self):
        """ building SN header

        |  Method to build header file from eqdsk with one single null (two points in PFxx) 
        |  The first five values of the eqdsk (nShot, tShot, modflg, FPPkat, IpiFPP) are already set correctly  
        |  The last quantities (rhoPF, PFL, Vol, Area, Qpl) are already set
        
        Parameters:
            None
        Attributes:
            None

        """

        print("Build hdr (SN)")

        nrho = len(self.eqdsk.rhopsi)
        dummy=np.linspace(0,1,nrho)
        
        self.hdr={'nSHOT':0,'tSHOT':0,'modflg':0,'FPPkat':0,'IpiFPP':self.eqdsk.Ip,\
                  'PFxx':[],'RPFx':[],'zPFx':[],'SSQ':[], 'devnam':self.devnam,\
                  'rhoPF':nrho,'PFL':dummy,'Vol':dummy,'Area':dummy,'Qpl':dummy} 

        #Find x-point
        f = plt.figure()
        ax2d = f.add_subplot(111)
        r,z = self.R_eqd, self.Z_eqd
        ax2d.contour(r,z, self.eqdsk.psi, 50)
        ax2d.set_title('choose x point position')
        x0 = plt.ginput()
        plt.close(f)
        self.xpoint = self._min_grad(x0=x0)        
        self.xflux = self.psi_coeff(self.xpoint[0], self.xpoint[1])*(2*math.pi)
        
        # find axis
        self.ax = self._min_grad(x0=[self.eqdsk.Raxis, self.eqdsk.Zaxis])     
        self.axflux = self.psi_coeff(self.ax[0], self.ax[1])*(2*math.pi)
        print("remember: I am multiplying psi axis and x-point times 2pi since in ascot it divides by it!")

        # poloidal flux of the special points.
        self.hdr['PFxx'] = [self.axflux[0][0], self.xflux[0][0]]
        self.hdr['RPFx'] = [self.ax[0], self.xpoint[0]]
        self.hdr['zPFx'] = [self.ax[1], self.xpoint[1]]
        self.hdr['SSQ']  = [self.eqdsk.R0EXP, self.eqdsk.Zaxis, 0, 0]
        
    def build_bkg(self):
        """ build bkg

        |  Method to build background file from eqdsk 
        |  -The first five values of the eqdsk (nShot, tShot, modflg, FPPkat, IpiFPP) are already set correctly  
        |  -The last quantities (rhoPF, PFL, Vol, Area, Qpl) are already set

        Parameters:
            None
        Attributes:
            None        
        
        """
        try:
            self.Fgrid.mean()
            print("Bphi already built!")
        except:
            self.calc_field()

        print("Build bkg")

        R_temp = np.linspace(self.eqdsk.rboxleft, self.eqdsk.rboxleft+self.eqdsk.rboxlength, self.nR)
        z_temp = np.linspace(-self.eqdsk.zboxlength/2., self.eqdsk.zboxlength/2., self.nz)
        #R_temp = np.linspace(float(np.around(np.min(self.R_w), decimals=2)), float(np.around(np.max(self.R_w), decimals=2)), self.nR)
        #z_temp = np.linspace(float(np.around(np.min(self.z_w), decimals=2)), float(np.around(np.max(self.z_w), decimals=2)), self.nz)

        psitemp = self.psi_coeff(z_temp, R_temp)
        #psitemp = self.psi_coeff(R_temp, z_temp)

        bphitemp = self.param_bphi(R_temp, z_temp)

        self.bkg={'type':'magn_bkg', 'phi0':0, 'nsector':0, 'nphi_per_sector':1,\
                  'ncoil':0, 'zero_at_coil':1,\
                  'R':R_temp,'z':z_temp, \
                  'phimap_toroidal':0, 'phimap_poloidal':0, \
                  'psi':[],\
                  'Bphi':bphitemp, 'BR':self.Br, 'Bz':self.Bz} 

        self.bkg['psi'] = psitemp*2*math.pi #in ASCOT Bfield, the psi is divided by 2*pi and reverses sign. This prevents it from happening  
        print("remember: I am multiplying psi times 2pi since in ascot it divides by it!")
    
    def _calc_psi_deriv(self):
        """ derivative of psi

        Compute the derivative of psi(poloidal flux) on a refined grid which will be used then for computing of the radial and vertical component of the magnetic field. It can be done by computing on a finer grid (128x128) within the separatrix

        Parameters:
            None
        Attributes:
            None
        """
        psi = self.eqdsk.psi
        self.dpsidR = np.zeros((self.eqdsk.nzbox, self.eqdsk.nrbox))
        self.dpsidZ = np.zeros((self.eqdsk.nzbox, self.eqdsk.nrbox))
        
        deriv = np.gradient(psi)
        # Note np.gradient gives y
        # derivative first, then x derivative
        ddR = deriv[1]
        ddZ = deriv[0]
        dRdi = np.asarray(1.0)/np.gradient(self.R_eqd)
        dRdi = np.tile(dRdi, [self.eqdsk.nzbox,1])
        dZdi = np.asarray(1.0)/np.gradient(self.Z_eqd)
        dZdi = np.tile(dZdi, [self.eqdsk.nrbox,1])
        dZdi = np.transpose(dZdi)
        #print("shape ddR:",np.shape(ddR),'shape dRdi:', np.shape(dRdi))
        #print('shape ddZ:',np.shape(ddZ),'shape dZdi:', np.shape(dZdi))
    
        self.dpsidR[:, :] = ddR*dRdi
        self.dpsidZ[:, :] = ddZ*dZdi


    def _min_grad(self, x0):
        """ minimum gradient

        find the point where there is the minimum of the flux

        Parameters:
            x0 (array): x,z coordinates of the starting point
        Attributes:
            None

        """
        try:
            self.dpsidR.mean()
        except:
            self._calc_psi_deriv()

        sp_dr = self.dpsidR
        sp_dz = self.dpsidZ
        R = self.R_eqd
        z = self.Z_eqd

        val_dr = interp.interp2d(R, z, sp_dr)
        val_dz = interp.interp2d(R, z, sp_dz)
        fun= lambda x: val_dr(x[0], x[1])**2 +val_dz(x[0], x[1])**2
        x = scipy.optimize.fmin(fun, x0)
        R0 = x[0]
        z0 = x[1]
        return R0,z0

    def calc_field(self):
        """ calculating magnetic field
        Function to calculate toroidal field (fields on poloidal plane set to 0) 

        Parameters:
            None
        Attributes:
            None

        """

        print("Calculating Bphi")
        inv_R = np.asarray(1.0)/np.array(self.R_eqd)
        inv_R = np.tile(inv_R,[self.eqdsk.nzbox, 1])
        self.Bphi = np.zeros(np.shape(inv_R))
        #Bphi is used, BR and Bz not but you must initialise it to 0 and
        # print them anyway
        #self.Br_t = -self.dpsidZ*inv_R
        #self.Bz_t =  self.dpsidR*inv_R
        self.Br = np.zeros((self.nR, self.nz))
        self.Bz = np.zeros((self.nR, self.nz))
        
        #Creating rhogrid and then Fgrid
        psinorm_grid = (self.eqdsk.psi-self.eqdsk.psiaxis)/(self.eqdsk.psiedge-self.eqdsk.psiaxis)
        self.rhogrid = np.sqrt(psinorm_grid)
        
        Fgrid=griddata(self.eqdsk.rhopsi, self.eqdsk.T, self.rhogrid, method='nearest')
        #Set values out of the separatrix (where Fgrid is NaN) to the value at the separatrix
        Fgrid[np.where(self.rhogrid>1.)] = self.eqdsk.T[-1]

        self.Fgrid=Fgrid 
        Bphi = np.multiply(Fgrid,inv_R)        
        self.param_bphi = interp.interp2d(self.R_eqd, self.Z_eqd, Bphi)
        self.Bphi = Bphi

    def write_head(self):
        """ writing header
        Write to input.magn_header file
        
        Parameters:
            None
        Attributes:
            None

        """
        try:
            hdr=self.hdr
        except:
            print("Build header first!")
            raise ValueError

        out_fname = 'input.magn_header'
        if self.devnam=='TCV':
            out_fname += '_'+self.infile[6:18]
			
        print('OUT header '+out_fname)
        outfile = open(out_fname, 'wa')
       
        
        #outfile.write('{:d} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n'.format(len(lines)))
        # shot info
        outfile.write('{:8d} {:10f} {:2d}\n'.format(hdr['nSHOT'], hdr['tSHOT'], hdr['modflg']))
        #device name        
        outfile.write(hdr['devnam'] +'\n')
        # something + plasma current        
        outfile.write('{:4d}   {:10f}\n'.format(hdr['FPPkat'], hdr['IpiFPP']))
        outfile.write('{:4d}\n'.format(len(hdr['PFxx'])))
        # Write the special points
        for j in range(len(hdr['PFxx'])):
            # poloidal flux
            outfile.write('{:8.6f} '.format(hdr['PFxx'][j]))
        outfile.write(' \n')

        for j in range(len(hdr['PFxx'])):
            # R
            outfile.write('{:8.6f} '.format(hdr['RPFx'][j]))
        outfile.write(' \n')
  
        for j in range(len(hdr['PFxx'])):
            # z
            outfile.write('{:8.6f} '.format(hdr['zPFx'][j]))
        outfile.write(' \n')
        
        #SSQ
        for i in xrange(0,len(hdr['SSQ']),4):
            tmp_str = ['{:8.6f} '.format(j) for j in hdr['SSQ'][i:i+4]]
            outfile.write(" ".join(tmp_str))
            outfile.write("\n")
        
        #print rhoPF 
        outfile.write(str(hdr['rhoPF'])+'\n')
        # other arrays
        
        for arr_name in ('PFL','Vol','Area','Qpl'):
            print("Writing ", arr_name)
            arr = hdr[arr_name]
            for i in xrange(0,len(arr),4):
                tmp_str = ['{:18.10f}'.format(j) for j in arr[i:i+4]]
                outfile.write(" ".join(tmp_str))
                outfile.write("\n")
        outfile.close()
        

    def write_bkg(self):
        """ write bkg
        Write to input.magn_bkg file
        
        self.bkg={'type':'magn_bkg', 'phi0':0, 'nsector':0, 'nphi_per_sector':1,\
                  'ncoil':18, 'zero_at_coil':1,\
                  'R':self.eqdsk.R,'z':self.eqdsk.Z, \
                  'phimap_toroidal':0, 'phimap_poloidal':0, \
                  'psi':-2*3.14*self.eqdsk.psi,\
                  'Bphi':self.Bphi, 'BR':self.Br, 'Bz':self.Bz}
       
        Parameters:
            None
        Attributes:
            None

        """
        try:
            self.bkg['Bphi'].mean()
        except:
            self.build_bkg()
            
        bkg=self.bkg
        out_fname = 'input.magn_bkg'
        if self.devnam=='TCV':
            out_fname += '_'+self.infile[6:18]

        print('OUT bkg '+out_fname)
        outfile = open(out_fname, 'wa') 
    
        #outfile.write('{:d} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n'.format(len(lines)))        
        outfile.write('{:18.10f} {:3d} {:3d} {:3d} {:3d}\n'.format(\
            bkg['phi0'], bkg['nsector'], bkg['nphi_per_sector'], bkg['ncoil'], \
            bkg['zero_at_coil']))
            
        outfile.write('{:18.10f} {:18.10f} {:3d}\n'.format(\
            bkg['R'][0], bkg['R'][-1], len(bkg['R'])))
        outfile.write('{:18.10f} {:18.10f} {:3d}\n'.format(\
            bkg['z'][0], bkg['z'][-1], len(bkg['z'])))
            
        if bkg['nsector'] ==0:
            # Do domething else if it's a different case,
            # but i didn't fully understand, anyway it's not my case yet
            outfile.write('{:d}\n'.format(0))
            outfile.write('{:d}\n'.format(0))
        else:
            print("Bkg[nsector] = ", bkg['nsector'])
            for arr_name in ('phimap_toroidal', 'phimap_poloidal'):
                arr = bkg[arr_name]
                for i in xrange(0,len(arr),18):
                    tmp_str = ['{:d}'.format(j) for j in arr[i:i+18]]
                    outfile.write(" ".join(tmp_str))
                    outfile.write("\n")
            
           
            #Bphi is used, BR and Bz not but you must initialise it to 0 and
           # print them anyway           
        for arr_name in ('psi', 'BR', 'Bphi', 'Bz'):
            print("Writing ", arr_name)
            arr_t = bkg[arr_name]
            #arr = self._perm_dims(arr_t)
            arr=arr_t

            #making the array plain:
            arr = arr.reshape(arr.size)

            for i in xrange(0,np.size(arr)-np.mod(np.size(arr),4),4):
                tmp_str = ['{:18.10f} {:18.10f} {:18.10f} {:18.10f}'.format(arr[i],arr[i+1],arr[i+2],arr[i+3])]
                outfile.write(" ".join(tmp_str))
                outfile.write("\n")
            tmp_str = ''
            for j in arr[-np.mod(np.size(arr),4):]:
                tmp_str += '{:18.10f} '.format(j)
            outfile.write(tmp_str)
            outfile.write("\n")           
                
        # Missing perturbation field, up to now useless
        
        outfile.close()
                 
    def _perm_dims(self,arr):
        """ permutating the dimensions
        This permutation of the array has to be done to correctly feed the input to ascot
        
        Parameters:
            arr (array): input array to permute
        Attributes:
            None

        """
        out_arr = []
        if len(np.shape(arr))==2:
            out_arr = np.transpose(arr)
        if len(np.shape(arr))==3:
            out_arr = np.transpose(arr, (1, 2, 0))
                
        return out_arr
            

    def _read_wall(self):
        """ read wall
        Reads 2D (R,Z) wall depending on the device name
        
        Parameters:
            None
        Attributes:
            None
        """
        try:
            if self.devnam == 'JT60SA':
                fname = '/home/vallar/JT60-SA/PARETI_2D_SA/input.wall_2d'
                #fname = "/home/vallar/JT60-SA/PARETI_2D_SA/input.wall_2d_clamped"
            elif self.devnam == 'TCV':
                fname = '/home/vallar/TCV/input.wall_2d'
                #fname = '/home/vallar/TCV/TCV_vessel_coord.dat'
            wall = np.loadtxt(fname, skiprows=1)
            self.R_w = wall[:,0]
            self.z_w = wall[:,1]
        except:
            print("No wall to read")
            self.R_w=[0]
            self.z_w=[0]
            return

class venus_Bfield(Bfield_eqdsk):
    """
    For the input given from L. Stipani, e.g. 53454 from VENUS/LEVIS
    """
    def __init__(self, infile_name, nR,nz):
        Bfield_eqdsk.__init__(self, infile_name, nR, nz)
        infile = sio.loadmat(infile_name)

        eq = infile['equilibrium']
        R_psi = eq['R'][0,0]#this is a 2D array with (s=rhotor**2, R)
        z_psi = eq['Z'][0,0]#this is a 2D array with (s=rhotor**2, Z)
        s_arr = eq['s'][0,0][0,:]
        r_arrt = eq['rplas'][0,0][0,:]
        r_arr = np.linspace(np.min(r_arrt), np.max(r_arrt), 30)
        z_arrt = eq['zplas'][0,0][0,:]
        z_arr = np.linspace(np.min(z_arrt), np.max(z_arrt), 30)

        psirz = np.zeros((len(r_arr), len(z_arr)), dtype=float)
        for i_s in np.arange(len(s_arr), step=1, dtype=int):
            R = R_psi[i_s,:]
            z = z_psi[i_s,:]
            for i_r in range(len(R)):
                indR = np.argmin(np.abs(r_arr-R[i_r]))
                indz = np.argmin(np.abs(z_arr-z[i_r]))
                if i_s==50000:
                    plot(r_arr[indR], z_arr[indz], 'x')
                psirz[indR, indz] = s_arr[i_s]
        psirz = np.transpose(psirz)

        psirz2 = psirz
        for ir,rttt in enumerate(r_arr):
            RLCFS_pos = r_arrt[z_arrt>0]
            indLCFS=np.argmin(np.abs(rttt-RLCFS_pos))
            indz = z_arr>z_arrt[indLCFS]
            psirz2[ir, indz] = 0
            
            RLCFS_neg = r_arrt[z_arrt<=0]
            indLCFS=np.argmin(np.abs(rttt-RLCFS_neg))
            indz = z_arr>z_arrt[indLCFS]
            psirz2[ir, indz] = 0
            
        params_psi = interp.interp2d(r_arr, z_arr, psirz, kind='linear')
        r_new = np.linspace(np.min(r_arrt), np.max(r_arrt), 200)
        z_new = np.linspace(np.min(z_arrt), np.max(z_arrt), 200)
        plt.contourf(r_new, z_new,params_psi(r_new, z_new), 20)
        
