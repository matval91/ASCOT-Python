"""
matteo.vallar@igi.cnr.it
Classes to handle magnetic fields I/O for ascot&bbnbi

"""
from __future__ import print_function

import numpy as np
import h5py, math
import matplotlib.pyplot as plt
import ReadEQDSK
from scipy.interpolate import griddata
import scipy.optimize
import scipy.interpolate as interp
import a4py.utils.cocos_transform as cocos
import a4py.plot.input.plot_input as plot_input

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
                        "q":"/boozer/qprof", "psi_rho":"/boozer/psi",\
                        "psi_sepa":"/boozer/psiSepa"}
        
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
            self.vardict[k] = self.infile[self.labdict[k]].value

        self.nrho = self.vardict['psi_rho'].shape[0]
        
        self.R = self.vardict['r'][:]
        self.z = self.vardict['z'][:]
        self.psi = self.vardict['psi_2D']
        self.psiedge = self.vardict['psi_sepa']
        self.rhopsi = self.vardict['psi_rho'][:]**0.5
        self.q = self.vardict['q']

        self.Bphi = self.vardict['bphi']

        
        #these nR and nZ are from the group /bfield/, the ones for /boozer/ are on rho and theta (spacing of grid)
        self.nR = self.vardict['nr'][:]        
        self.nZ = self.vardict['nz'][:]

        self.rmin=np.min(self.R)
        self.rmax=np.max(self.R)
        self.zmin=np.min(self.z)
        self.zmax=np.max(self.z)

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
           
            
    def plot_B(self,f=0 ):
        plot_input.plot_Bfield(self,f)
        
        
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
        self.R = self.R_eqd
        self.z = self.Z_eqd
        self.psi = self.eqdsk.psi
        self.psiedge = self.eqdsk.psiedge
        self.rhopsi = self.eqdsk.rhopsi
        self.q = self.eqdsk.q
        self.nR=len(self.R)
        self.nz=len(self.z)
        self.nR=nR
        self.nz=nz        

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
        self._cocos_transform(self.COCOS)
        #This is for Equilibrium from CREATE for scenario 2, also to modify in build bkg
        self.psi_coeff = interp.RectBivariateSpline(self.Z_eqd, self.R_eqd, self.eqdsk.psi)   
        
    def _cocos_transform(self, COCOS):
        """ cocos transformations
        This function converts the magnetic input from their starting cocos to eqdsk 5 (needed by ascot).

        Parameters:
            COCOS (int): input cocos. Now useable only 2,3,4,7,12,13,14,17
        Attributes:
            None
        
        """
        cocos.cocos_transform(self.eqdsk, COCOS, 5, \
                              sigma_ip_out=1, sigma_b0_out=-1)

    def plot_B(self, f=0):
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
            
        plot_input.plot_Bfield(self,f)

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
        ax2d.axis('equal')
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
        self.hdr['RPFx'] = [self.xpoint[0], self.ax[0]]
        self.hdr['zPFx'] = [self.xpoint[1], self.ax[1]]
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
            self.param_bphi.x
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
                  'Bphi':bphitemp, 'BR':self.Br, 'Bz':self.Bz, \
                  'Bphi_pert':self.Bphi_pert, 'BR_pert':self.BR_pert, 'Bz_pert':self.Bz_pert} 

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

        self._Fgrid=Fgrid 
        Bphi = np.multiply(Fgrid,inv_R)        
        self.param_bphi = interp.interp2d(self.R_eqd, self.Z_eqd, Bphi)
        self.Bphi = Bphi
        # Adding perturbation fields
        self.Bphi_pert = self.Br
        self.BR_pert = self.Br
        self.Bz_pert = self.Br
        
        
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
        outfile = open(out_fname, 'w')
       
        
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
        for i in range(0,len(hdr['SSQ']),4):
            tmp_str = ['{:8.6f} '.format(j) for j in hdr['SSQ'][i:i+4]]
            outfile.write(" ".join(tmp_str))
            outfile.write("\n")
        
        #print rhoPF 
        outfile.write(str(hdr['rhoPF'])+'\n')
        # other arrays
        
        for arr_name in ('PFL','Vol','Area','Qpl'):
            print("Writing ", arr_name)
            arr = hdr[arr_name]
            for i in range(0,len(arr),4):
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
        outfile = open(out_fname, 'w') 
    
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
                for i in range(0,len(arr),18):
                    tmp_str = ['{:d}'.format(j) for j in arr[i:i+18]]
                    outfile.write(" ".join(tmp_str))
                    outfile.write("\n")
            
           
            #Bphi is used, BR and Bz not but you must initialise it to 0 and
           # print them anyway           
        for arr_name in ('psi', 'BR', 'Bphi', 'Bz'): #, 'Bphi_pert', 'BR_pert', 'Bz_pert'):
            print("Writing ", arr_name)
            arr_t = bkg[arr_name]
            arr = arr_t;
            #if arr_name!='psi':arr = self._perm_dims(arr_t)

            #making the array plain:
            arr = arr.reshape(arr.size)
            for i in range(0,np.size(arr)-np.mod(np.size(arr),4),4):
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
