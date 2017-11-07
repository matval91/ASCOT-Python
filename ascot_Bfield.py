"""
Class for magnetic field
"""
import numpy as np
import h5py, math
import matplotlib.pyplot as plt
import ascot_misc

import ReadEQDSK_MV
import scipy.optimize
import scipy.interpolate as interp
#from matplotlib.figure import Figure
import scipy.io as sio

class Bfield_out:
    """
    Class for handling the magnetic field specifications and plots
    Groups in h5 file (09/01/2017, ascot4)
    /bfield                  Group
    /bfield/2d               Group
    /bfield/2d/bphi          Dataset {600, 300}
    /bfield/2d/psi           Dataset {600, 300}
    /bfield/nphi             Dataset {1}
    /bfield/nr               Dataset {1}
    /bfield/nz               Dataset {1}
    /bfield/r                Dataset {300}
    /bfield/raxis            Dataset {1}
    /bfield/z                Dataset {600}
    /bfield/zaxis            Dataset {1}


    /boozer                  Group
    /boozer/Babs             Dataset {200, 100}
    /boozer/Ifunc            Dataset {100}
    /boozer/R                Dataset {200, 100}
    /boozer/axisR            Dataset {1}
    /boozer/axisz            Dataset {1}
    /boozer/btheta           Dataset {200, 100}
    /boozer/delta            Dataset {200, 100}
    /boozer/gfunc            Dataset {100}
    /boozer/nu               Dataset {200, 100}
    /boozer/psi              Dataset {100}
    /boozer/psiAxis          Dataset {1}
    /boozer/psiMax           Dataset {1}
    /boozer/psiMin           Dataset {1}
    /boozer/psiSepa          Dataset {1}
    /boozer/qprof            Dataset {100}
    /boozer/rgridmax         Dataset {1}
    /boozer/rgridmin         Dataset {1}
    /boozer/theta            Dataset {200}
    /boozer/z                Dataset {200, 100}
    /boozer/zgridmax         Dataset {1}
    /boozer/zgridmin         Dataset {1}

    
    METHODS:
    __init__(self, infile) to store variables
    checkplot(self) to plot magnetic values and check them

    """



    
    def __init__(self, infile_name):
        #defining dictionary for handling data from h5file
        self.labdict = {"bphi":"/bfield/2d/bphi","psi_2D":"/bfield/2d/psi",\
                        "nr":"/bfield/nr", "nz":"/bfield/nz",\
                        "r":"/bfield/r", "z":"/bfield/z",\
                        "raxis":"/bfield/raxis", "zaxis":"/bfield/zaxis",\
						\
                        "Babs":"/boozer/Babs", "I":"/boozer/Ifunc",\
                        "r_boo":"/boozer/R","axisR_boo":"/boozer/axisR", "axisZ_boo":"/boozer/axisz",\
                        "btheta":"/boozer/btheta", "delta":"/boozer/delta",\
                        "g":"/boozer/gfunc", "nu":"/boozer/nu",\
                        "psi_rho":"/boozer/psi", "psi_ax":"/boozer/psiAxis",\
                        "maxPsi":"/boozer/psiMax", "minPsi":"/boozer/psiMin",\
                        "psi_sepa":"/boozer/psiSepa","q":"/boozer/qprof",\
                        "rmax":"/boozer/rgridmax", "rmin":"/boozer/rgridmin",\
                        "theta":"/boozer/theta",\
                        "z_boo":"/boozer/z","zmax":"/boozer/zgridmax", "zmin":"/boozer/zgridmin"\
                        }

        infile=h5py.File(infile_name)

        self.vardict = {}
        #store data with the correct labels
        
        for k in self.labdict.keys():
            self.vardict[k] = infile[self.labdict[k]]

        self.nrho = self.vardict['psi_rho'].shape[0]
        self.rho = np.zeros(self.nrho)
        self.rho = self.vardict['psi_rho'][:]**0.5
        self.nthe = self.vardict['theta'].shape[0]
        self.the = np.zeros(self.nthe)
        self.the = self.vardict['theta'][:]


##         self.rmin = self.vardict["rmin"][:]
##         self.zmin = self.vardict["zmin"][:]
##         self.rmax = self.vardict["rmax"][:]
##         self.zmax = self.vardict["zmax"][:]

        #these nR and nZ are from the group /bfield/, the ones for /boozer/ are on rho and theta (spacing of grid)
        self.nR = self.vardict['nr'][:]        
        self.nZ = self.vardict['nz'][:]
        self.R = self.vardict['r'][:]
        self.Z = self.vardict['z'][:]
        self.rmin=np.min(self.R)
        self.rmax=np.max(self.R)
        self.zmin=np.min(self.Z)
        self.zmax=np.max(self.Z)
        #self.R = np.linspace(self.rmin, self.rmax, self.nR)
        #self.Z = np.linspace(self.zmin, self.zmax, self.nZ)

        self.w = ascot_misc.misc(infile_name)
        
        # doing plot to check what we have
        #self.checkplot()


    def checkplot(self):
        """
        Method to plot the values and check the magnetic field we are looking at
        
        """

        n_row = 1
        n_col = 3
        for e,val in enumerate(["psi_2D","q","bphi"]):

            plt.subplot(n_row,n_col,e+1)
            plt.title(val)

            #plot for mag surfaces 2D
            if val == "psi_2D":
                r,z = self.vardict['r'], self.vardict['z']
                plt.contour(r,z,self.vardict[val], self.nrho)
                plt.contour(r,z,self.vardict[val],1)
                #CB = plt.colorbar(CS, shrink=0.8, extend='both')
                plt.xlabel("R [m]")
                plt.ylabel("Z [m]")
                plt.plot(self.w.vardict["R_wall"], self.w.vardict["Z_wall"], 'k')

            elif val == 'bphi': #plot of g
                plt.plot(self.vardict['r'][:], self.vardict[val][0,:]*self.vardict['r'][:])
                plt.xlabel(r'R[m]')
                plt.ylabel(r'Rx$B_\phi$')
            elif val == 'q': #plot of q
                plt.plot(self.rho, self.vardict[val])
                plt.xlabel(r'$\rho_{POL}$')
                plt.ylabel(r'q')
        plt.show()


        
class Bfield_in:
    """
    Script for writing and reading the magnetic background
    porting from matlab (16-1-2017)
    
    """
    
    def __init__(self, infile, nR, nz):
        """
        Initialisation, with EQDSK:
        -EQDSK (infile = eqdsk file)
        """
        self.import_from_eqdsk(infile)
        #these are the dimensions of the output arrays
        self.nR=nR
        self.nz=nz
        

    def import_from_eqdsk(self, infile_eqdsk):
        """
        function for import from eqdsk file
        this is the structure of the eqdsk struct:
        
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
    
    
    
        self.eqdsk= ReadEQDSK_MV.ReadEQDSK(infile_eqdsk)
        self.eqdsk.psi = np.reshape(self.eqdsk.psi, (self.eqdsk.nzbox, self.eqdsk.nrbox))
        #self.eqdsk.psi = -2.*math.pi*self.eqdsk.psi
        self.psi_coeff = interp.interp2d(self.eqdsk.R_grid, self.eqdsk.Z_grid, self.eqdsk.psi)
        self._read_wall()


    def eqdsk_checkplot(self):
        """
        Method to plot the values and check the magnetic field we are looking at
        
        """
        n_row = 1
        n_col = 3

        f = plt.figure()
        ax2d = f.add_subplot(131)
        r,z = self.eqdsk.R_grid, self.eqdsk.Z_grid
        ax2d.contour(r,z, self.eqdsk.psi, 50)
        ax2d.contour(r,z, self.eqdsk.psi, [1], linewidth=3, linecolor='k')
        #CB = plt.colorbar(CS, shrink=0.8, extend='both')
        ax2d.set_xlabel("R")
        ax2d.set_ylabel("Z")
        if self.R_w[0]!=0:
            ax2d.plot(self.R_w, self.z_w, 'k',linewidth=2)

        axq = f.add_subplot(132)
        axq.plot(self.eqdsk.rhopsi, self.eqdsk.q)
        axq.set_xlabel(r'$\rho_{POL}$')
        axq.set_ylabel(r'q')

        axf = f.add_subplot(133)
        axf.plot(self.eqdsk.R_grid, self.eqdsk.T)
        axf.set_xlabel(r'R [m]')
        axf.set_ylabel(r'g (poloidal flux)')

        plt.show()


    def write(self):
        self.write_head()
        self.write_bkg()

    def build_lim(self):
        self.build_header_lim()
        self.build_bkg()    

    def build_SN(self):
        self.build_header_SN()
        self.build_bkg()
        
    def build_header_lim(self):
        """
        Method to build header file from eqdsk without single nulls (one point in PFxx) 
        -The first five values of the eqdsk (nShot, tShot, modflg, FPPkat, IpiFPP)
        are already set correctly  
        -The last quantities (rhoPF, PFL, Vol, Area, Qpl) are already set
        
        """
        print "Build hdr (limiter)"

        
        nrho = len(self.eqdsk.rhopsi)
        dummy=np.linspace(0,1,nrho)
        
        self.hdr={'nSHOT':0,'tSHOT':0,'modflg':0,'FPPkat':0,'IpiFPP':self.eqdsk.Ip,\
                  'PFxx':[],'RPFx':[],'zPFx':[],'SSQ':[], 'devnam':'JT-60SA',\
                  'rhoPF':129,'PFL':dummy,'Vol':dummy,'Area':dummy,'Qpl':dummy} 
        #find derivatives for finding x point
        self.dpsidR, self.dpsidZ = self._calc_psi_deriv()  
        
        # find axis
        self.ax = self.min_grad(x0=[self.eqdsk.Raxis, self.eqdsk.Zaxis])     
        self.axflux = self.psi_coeff(self.ax[0], self.ax[1])*(2*math.pi)

        # poloidal flux of the special points (only one in this case)
        self.hdr['PFxx'] = [self.axflux[0]]#-self.axflux[0]-self.xflux[0]
        self.hdr['RPFx'] = [self.ax[0]]
        self.hdr['zPFx'] = [self.ax[1]]
        self.hdr['SSQ']  = [self.eqdsk.R0EXP, self.eqdsk.Zaxis, 0, 0]

        
    def build_header_SN(self):
        """
        Method to build header file from eqdsk with one single null (two points in PFxx) 
        -The first five values of the eqdsk (nShot, tShot, modflg, FPPkat, IpiFPP)
        are already set correctly  
        -The last quantities (rhoPF, PFL, Vol, Area, Qpl) are already set
        
        """

        print "Build hdr (SN)"

        nrho = len(self.eqdsk.rhopsi)
        dummy=np.linspace(0,1,nrho)
        
        self.hdr={'nSHOT':0,'tSHOT':0,'modflg':0,'FPPkat':0,'IpiFPP':self.eqdsk.Ip,\
                  'PFxx':[],'RPFx':[],'zPFx':[],'SSQ':[], 'devnam':'JT-60SA',\
                  'rhoPF':129,'PFL':dummy,'Vol':dummy,'Area':dummy,'Qpl':dummy} 
        #Smooth psi and find derivatives for finding x point
        #self.psi_sm=self.smooth_psi()
        self.dpsidR, self.dpsidZ = self._calc_psi_deriv()  

        #Find x-point
        f = plt.figure()
        ax2d = f.add_subplot(111)
        r,z = self.eqdsk.R_grid, self.eqdsk.Z_grid
        ax2d.contour(r,z, self.eqdsk.psi, 50)
        ax2d.set_title('choose x point position')
        x0 = plt.ginput()
        self.xpoint = self.min_grad(x0=x0)        
        self.xflux = self.psi_coeff(self.xpoint[0], self.xpoint[1])
        
        # find axis
        self.ax = self.min_grad(x0=[self.eqdsk.Raxis, self.eqdsk.Zaxis])     
        self.axflux = self.psi_coeff(self.ax[0], self.ax[1])

        # poloidal flux of the special points.
        self.hdr['PFxx'] = [self.xflux[0], self.axflux[0]]#-self.axflux[0]-self.xflux[0]
        self.hdr['RPFx'] = [self.xpoint[0], self.ax[0]]
        self.hdr['zPFx'] = [self.xpoint[1], self.ax[1]]
        self.hdr['SSQ']  = [self.eqdsk.R0EXP, 0, 0, 0]
        
        #THESE HAVE BEEN ALREADY INITIALISED
        #rubbish = np.linspace(0,1,250)
        #self.hdr['rhoPF'] = rubbish
        #self.hdr['PFL']   = rubbish
        #self.hdr['Vol']   = rubbish
        #self.hdr['Area']  = rubbish
        #self.hdr['Qpl']   = rubbish

    def build_bkg(self):
        """
        Method to build background file from eqdsk 
        -The first five values of the eqdsk (nShot, tShot, modflg, FPPkat, IpiFPP)
        are already set correctly  
        -The last quantities (rhoPF, PFL, Vol, Area, Qpl) are already set
        
        """
        try:
            self.param_bphi(self.eqdsk.R_grid, self.eqdsk.Z_grid)
        except:
            self.calc_field()

        print "Build bkg"

        R_temp=np.linspace(min(self.eqdsk.R_grid), max(self.eqdsk.R_grid),self.nR)
        z_temp=np.linspace(min(self.eqdsk.Z_grid), max(self.eqdsk.Z_grid),self.nz)
        psitemp = self.psi_coeff(R_temp, z_temp)
        bphitemp = self.param_bphi(R_temp, z_temp)
        self.bkg={'type':'magn_bkg', 'phi0':0, 'nsector':0, 'nphi_per_sector':1,\
                  'ncoil':0, 'zero_at_coil':1,\
                  'R':R_temp,'z':z_temp, \
                  'phimap_toroidal':0, 'phimap_poloidal':0, \
                  'psi':[],\
                  'Bphi':bphitemp, 'BR':self.Br, 'Bz':self.Bz} 

        self.bkg['psi'] = psitemp
        
#    def smooth_psi(self):
#        psi=self.eqdsk.psi
#        R=self.eqdsk.R_grid
#        z=self.eqdsk.Z_grid
#        bbox=[min(R), max(R), min(z), max(z)]
##        newpsi=interp.SmoothBivariateSpline(x=R, y=z, z=psi, bbox=bbox, kx=3, ky=3)
#        newpsi_coeff = interp.bisplrep(R,z,psi)
#       
#        return newpsi        
        
    def _calc_psi_deriv(self):
        """
        Compute the derivative of the poloidal flux on a refined grid which
        will be used then for computing of the radial and vertical component of
        the magnetic field.
        It can be done by computing on a finer grid (128x128)
        within the separatrix
        for each time
        """
        psi = self.eqdsk.psi
        dpsidR = np.zeros((self.eqdsk.nzbox, self.eqdsk.nrbox))
        dpsidZ = np.zeros((self.eqdsk.nzbox, self.eqdsk.nrbox))
        
        #dpsidR = np.zeros((np.size(self.eqdsk.R_grid), np.size(self.eqdsk.Z_grid)))
        #dpsidZ = np.zeros((np.size(self.eqdsk.R_grid), np.size(self.eqdsk.Z_grid)))
        deriv = np.gradient(psi)
        # Note np.gradient gives y
        # derivative first, then x derivative
        ddR = deriv[1]
        # ddR = self.psi(Rgrid,Zgrid,dx=1)
        ddZ = deriv[0]
        # ddZ = self.psi(Rgrid,Zgrid,dy=1)
        dRdi = 1.0/np.gradient(self.eqdsk.R_grid)
        dRdi = np.tile(dRdi, [self.eqdsk.nzbox,1])
        #dRdi = np.transpose(dRdi)
        dZdi = 1.0/np.gradient(self.eqdsk.Z_grid)
        dZdi = np.tile(dZdi, [self.eqdsk.nrbox,1])
        dZdi = np.transpose(dZdi)
        #print "shape ddR:",np.shape(ddR),'shape dRdi:', np.shape(dRdi)
        #print 'shape ddZ:',np.shape(ddZ),'shape dZdi:', np.shape(dZdi)
    
        dpsidR[:, :] = ddR*dRdi
        dpsidZ[:, :] = ddZ*dZdi
    
        return dpsidR, dpsidZ   
        
    def min_grad(self, x0):
        sp_dr = self.dpsidR
        sp_dz = self.dpsidZ
        R = self.eqdsk.R_grid
        z = self.eqdsk.Z_grid

        val_dr = interp.interp2d(R, z, sp_dr)
        val_dz = interp.interp2d(R, z, sp_dz)
        fun= lambda x: val_dr(x[0], x[1])**2 +val_dz(x[0], x[1])**2
#        fun = np.multiply(sp_dr, sp_dr)+np.multiply(sp_dz, sp_dz)
        x = scipy.optimize.fmin(fun, x0)
        R0 = x[0]
        z0 = x[1]
        return R0,z0

    def calc_field(self):
        """
        Function to calculate fields 
        """
        try:
            self.dpsidR.mean()
        except:
            self.dpsidR, self.dpsidZ = self._calc_psi_deriv()

        print "Calculating Bphi"
        #tmpR = np.linspace(np.min(self.eqdsk.R_limits), np.max(self.eqdsk.R_limits), np.size(self.eqdsk.R_grid))
        inv_R = 1.0/np.array(self.eqdsk.R_grid)
        inv_R = np.tile(inv_R,[self.eqdsk.nzbox, 1])
        #inv_R = np.transpose(inv_R)
        #Bphi is used, BR and Bz not but you must initialise it to 0 and
        # print them anyway
        #self.Br = -self.dpsidR*inv_R
        #self.Bz =  self.dpsidZ*inv_R
        self.Br = np.zeros((self.nR, self.nz))
        self.Bz = np.zeros((self.nR, self.nz))
     
        f = np.array(self.eqdsk.T)
        f = f[np.newaxis,:]
        Bphi = f*inv_R

        self.param_bphi = interp.interp2d(self.eqdsk.R_grid,self.eqdsk.Z_grid, Bphi)
        #R_temp=np.linspace(min(self.eqdsk.R_limits), max(self.eqdsk.R_limits), self.nR)
        #z_temp=np.linspace(min(self.eqdsk.Z_limits), max(self.eqdsk.Z_limits), self.nz)
        #self.Bphi = param_bphi(R_temp, z_temp)
    
    def write_head(self):
        """
        Write to input.magn_header file
        """
        try:
            hdr=self.hdr
        except:
            self.build_header()
            hdr = self.hdr
        out_fname = 'input.magn_header'
        outfile = open(out_fname, 'wa')
        
        #frmstr1 = '%8.6g %8.6g %8.6g %8.6g\n';
        #frmstr2 = '%18.10g %18.10g %18.10g %18.10g\n';        
        
        
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
            print "Writing ", arr_name
            arr = hdr[arr_name]
            for i in xrange(0,len(arr),4):
                tmp_str = ['{:18.10f}'.format(j) for j in arr[i:i+4]]
                outfile.write(" ".join(tmp_str))
                outfile.write("\n")
        outfile.close()
        

    def write_bkg(self):
        """
        Write to input.magn_bkg file
        
        self.bkg={'type':'magn_bkg', 'phi0':0, 'nsector':0, 'nphi_per_sector':1,\
                  'ncoil':18, 'zero_at_coil':1,\
                  'R':self.eqdsk.R,'z':self.eqdsk.Z, \
                  'phimap_toroidal':0, 'phimap_poloidal':0, \
                  'psi':-2*3.14*self.eqdsk.psi,\
                  'Bphi':self.Bphi, 'BR':self.Br, 'Bz':self.Bz}
       
        """
        try:
            self.bkg['Bphi'].mean()
        except:
            self.build_bkg()
            
        bkg=self.bkg
        out_fname = 'input.magn_bkg'
        outfile = open(out_fname, 'wa')
        
        #frmstr1 = '%8.6g %8.6g %8.6g %8.6g\n';
        #frmstr2 = '%18.10g %18.10g %18.10g %18.10g\n';        
        
        
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
            print "Bkg[nsector] = ", bkg['nsector']
            for arr_name in ('phimap_toroidal', 'phimap_poloidal'):
                arr = bkg[arr_name]
                for i in xrange(0,len(arr),18):
                    tmp_str = ['{:d}'.format(j) for j in arr[i:i+18]]
                    outfile.write(" ".join(tmp_str))
                    outfile.write("\n")
            
           
            #Bphi is used, BR and Bz not but you must initialise it to 0 and
           # print them anyway           
        for arr_name in ('psi', 'BR', 'Bphi', 'Bz'):
            print "Writing ", arr_name
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
        out_arr = []
        if len(np.shape(arr))==2:
            out_arr = np.transpose(arr)
        if len(np.shape(arr))==3:
            out_arr = np.transpose(arr, (1, 2, 0))
                
        return out_arr
            

    def _read_wall(self):
        
        try:
            fname = "/home/vallar/JT60-SA/PARETI_2D_SA/input.wall_2d"
            fname = "/home/vallar/JT60-SA/PARETI_2D_SA/input.wall_2d_clamped"
            fname = '/home/vallar/TCV/from_jari/input.wall_2d'           
            wall = np.loadtxt(fname, skiprows=1)
            self.R_w = wall[:,0]
            self.z_w = wall[:,1]

        except:
            print "No wall to read"
            self.R_w=[0]
            self.z_w=[0]
            return

class venus_Bfield(Bfield_in):
    """
    For the input given from L. Stipani, e.g. 53454 from VENUS/LEVIS
    """
    def __init__(self, infile_name, nR,nz):
        Bfield_in.__init__(self, infile_name, nR, nz)
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
        
