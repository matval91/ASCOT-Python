"""
Class for profiles
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy import interpolate
import scipy.io as sio
import time
import math

colours = ['k','b','r','c','g',\
           'k','b','r','c','g',\
           'k','b','r','c','g']

class profiles:
    """
    Superclass with the following methods:
    -__init__: initialises the profiles we need
    -plot
    -write
    - _spline
    - _extrapolate
    """
    def __init__(self):
        self.rho=[]
        self.te=[]
        self.ti=[]
        self.ne=[]
        self.ni1=[]
        self.ni2=[]
        self.ni3=[]
        self.vtor=[]
        self.zeff=[]

        self.nion=1
        self.Z=[]
        self.A=[]
        self.coll_mode=[]
      

    def write_input(self):
        """
        Method to write the input.plasma_1d file
        """
        
        out_fname = "input.plasma_1d"
        outfile = open(out_fname, 'wa')
        
        outfile.write('# Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation \n')
        outfile.write('# range must cover [0,1] of normalised poloidal rho. It can exceed 1. \n')
        outfile.write('# {:s} (first 3 lines are comment lines) \n'.format(time.strftime('%d%b%y')))
            
        outfile.write('{:d}\t{:1d}\t# Nrad,Nion \n'.format(self.nrho,self.nion))
        strcoll = str(1)+' ' # for electrons
        strZ=''
        strA=''
        for i in range(self.nion):
            strZ += str(self.Z[i]) + ' '
            strA += str(self.A[i]) + ' '
            strcoll += str(int(self.coll_mode[i])) + ' '
        strZ +='\t\t# ion Znum \n'
        strA +='\t\t# ion Amass \n'
        strcoll += '# collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons \n'
        outfile.write(strZ)				
        outfile.write(strA)
        outfile.write(strcoll)
            
        lab_len=15
        strlabel='RHO (pol)'.ljust(lab_len)+'Te (eV)'.ljust(lab_len)+'Ne (1/m3)'.ljust(lab_len)+'Vtor_I (rad/s)'.ljust(lab_len)+\
                  'Ti1 (eV)'.ljust(lab_len)
        for i in range(self.nion):
            tmpstr ='Ni{:d} (1/m3)'.format(i+1)
            strlabel+=tmpstr.ljust(lab_len)
        strlabel+='\n'
        outfile.write(strlabel)
        data=np.array((self.rho, self.te, self.ne, self.vt, self.ti), dtype=float)
        data = np.concatenate([data, [self.ni[i,:] for i in range(self.nion)]])

        data=np.transpose(data)
        np.savetxt(outfile, data, fmt='%15e')
        outfile.close()

    def plot_profiles(self):
        #=====================================================================================
        # SET TEXT FONT AND SIZE
        #=====================================================================================
        #plt.rc('font', family='serif', serif='Palatino')
        #plt.rc('text', usetex=True)
        plt.rc('linethick',)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=====================================================================================
        fig=plt.figure()
        axte = fig.add_subplot(321)
        axne = fig.add_subplot(322)
        axti = fig.add_subplot(323)
        axni = fig.add_subplot(324)
        axvt = fig.add_subplot(325)
        lw=2.
        axte.plot(self.rho, self.te*1e-3,'k', linewidth=lw)
        axne.plot(self.rho, self.ne,'k', linewidth=lw)
        axti.plot(self.rho, self.ti*1e-3,'k', linewidth=lw)
        axvt.plot(self.rho, self.vt,'k', linewidth=lw)
        for i in range(self.nion):
            axni.plot(self.rho, self.ni[i,:],colours[i], label='A='+str(self.A[i]), linewidth=lw)
        
        #axni.legend(loc='best')

        for ax in [axte, axne, axti, axni, axvt]:
            ax.set_xlabel(r'$\rho$')
            #=====================================================================================
            # Create your ticker object with M ticks
            M = 4
            yticks = ticker.MaxNLocator(M)
            xticks = ticker.MaxNLocator(M)
            # Set the yaxis major locator using your ticker object. You can also choose the minor
            # tick positions with set_minor_locator.
            ax.yaxis.set_major_locator(yticks)
            ax.xaxis.set_major_locator(xticks)
            #=====================================================================================
            
        axte.set_ylabel(r'$T_e$ [keV]')
        axne.set_ylabel(r'$N_e$ [1/$m^3$]')
        axti.set_ylabel(r'$T_i$ [keV]')
        axni.set_ylabel(r'$N_i$ [1/$m^3$]')
        axvt.set_ylabel(r'$V_{tor} [rad/s]$')
        #fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
        plt.show()


    def _spline(self,rho,data, rho_new):
        """
        Function to evaluate splines of the data input, in order to make an output in the right rho grid
        """
        dummy    = interpolate.InterpolatedUnivariateSpline(rho, data, ext=0)
        data_new = dummy(rho_new)
        return data_new

    def _extrapolate(self):
        """
        Function that does an exponential fit over rho=1 surface. Is better to give it as input instead of letting it to ASCOT
        """
        x = np.linspace(1.001, 1.2, 200)
        dec_l = 0.01
        ni_ov = np.zeros((self.nion, len(x)), dtype=float)
        ninew = np.zeros((self.nion, self.nrho+len(x)),dtype=float)
        ne_ov1 = self.ne[self.nrho-2]*np.exp(-((x-1.)/dec_l))
        te_ov1 = self.te[self.nrho-2]*np.exp(-(x-1.)/dec_l)
        ti_ov1 = self.ti[self.nrho-2]*np.exp(-(x-1.)/dec_l)
        vt_ov1 = self.vt[self.nrho-2]*np.exp(-(x-1.)/dec_l)
        for i in range(self.nion):
            ni_ov[i,:] = self.ni[i,self.nrho-2]*np.exp(-(x-1.)/dec_l)
            ninew[i,:] = np.concatenate([self.ni[i,:], ni_ov[i,:]])
        self.ni = ninew
        self.rho = np.concatenate([self.rho, x])
        self.nrho += len(x)
        self.ne  = np.concatenate([self.ne, ne_ov1])
        self.te  = np.concatenate([self.te, te_ov1])
        self.ti  = np.concatenate([self.ti, ti_ov1])
        self.vt  = np.concatenate([self.vt, vt_ov1])

class h5_profiles(profiles):   
    """
    Class for handling the profiles data from h5 file
    
    DATA in h5 file (09/01/2017, ascot4)

    /plasma                  Group
    /plasma/1d               Group
    /plasma/1d/ne            Dataset {1373}
    /plasma/1d/ni            Dataset {1373, 3}
    /plasma/1d/rho           Dataset {1373}
    /plasma/1d/te            Dataset {1373}
    /plasma/1d/ti            Dataset {1373}
    /plasma/1d/vtor          Dataset {1373}
    /plasma/1d/zeff          Dataset {1373}
    /plasma/anum             Dataset {3}
    /plasma/colls            Dataset {4}
    /plasma/znum             Dataset {3}

    
    METHODS:
    __init__(self, infile) to store variables
    checkplot(self) to plot values and check them

    """

    def __init__(self, infile_name, **kwargs):
        profiles.__init__(self)
        #defining dictionary for handling data from h5file
        self.labdict = {"rho":"/plasma/1d/rho",\
                        "ne":"/plasma/1d/ne","ni":"/plasma/1d/ni",\
                        "te":"/plasma/1d/te","ti":"/plasma/1d/ti",\
                        "vtor":"/plasma/1d/vtor","zeff":"/plasma/1d/zeff",\
                        "a_ions":"/plasma/anum", "znum":"/plasma/znum"\
                        }
        self.inf_name = infile_name
        if infile_name[-2:]=='h5':
            self.read_h5()
           
            
    def read_h5(self):
        infile = h5py.File(self.inf_name,'r')

        vardict = self.labdict
        #store data with the correct labels
        for k in infile['plasma/1d'].keys():
            vardict[k] = infile[self.labdict[k]]
        vardict['a_ions']=infile['/plasma/anum'].value
        vardict['znum']=infile['/plasma/znum'].value
        

        rho = vardict['rho']
        self.nrho = np.size(self.rho)
        if vardict['a_ions'][0]!='/':
            self.nspec = len(vardict['a_ions'])
        else:
            self.nspec = vardict['ni'].shape[1]
        print "Number of ions: ", self.nspec
        if len(vardict['a_ions'])!=len(vardict['znum']):
            print "ERROR! array of A and Z don't have the same length"

        self.A = vardict['a_ions']
        self.Z = vardict['znum']
        self.nion = self.nspec
        
        self.rho_in = rho
        self.te_in  = vardict['te']
        self.ne_in  = vardict['ne']  
        self.ti_in  = vardict['ti']
        self.ni1_in  = vardict['ni'][:,0]
        if self.nion==2:
            self.ni2_in  = vardict['ni'][:,1]
        elif self.nion==3:
            self.ni2_in  = vardict['ni'][:,1]
            self.ni3_in  = vardict['ni'][:,2]
        self.vt_in  = vardict['vtor']
        self.zeff_in  = vardict['zeff']

        self.smooth()
        # plot to check what we have
        #self.checkplot()

    def smooth(self):
        """
        smooth input data to data wanted
        """
        n_rhonew = 1000
        self.rho = np.linspace(0,1,n_rhonew)
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        self.ni1 = self._spline(self.rho_in, self.ni1_in, self.rho)
        if self.nion==2:
            self.ni2 = self._spline(self.rho_in, self.ni2_in, self.rho)
        elif self.nion==3:
            self.ni2 = self._spline(self.rho_in, self.ni2_in, self.rho)
            self.ni2 = self._spline(self.rho_in, self.ni3_in, self.rho)

        self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        self.zeff = self._spline(self.rho_in, self.zeff_in, self.rho)

class dat_profiles(profiles):
    """
    Function to write the profile file for ASCOT from an ascii file in the format (rho, quantity)
    
    INPUT:
    -names is a list containing the name of the files for Te, Ne, vtor_i, ti
    These files should have one column w the rho and one with the data
    -nion is the number of ions (=1 if not set)
    
    OUTPUT:
    input.plasma_1d
    
    # Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation
    # range must cover [0,1] of normalised poloidal rho. It can exceed 1. 
    # 18Jan08 for testing (first 3 lines are comment lines)
    111	3	# Nrad,Nion
    1 1 6		# ion Znum				
    1 2 12		# ion Amass
    1 1 1 1		# collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons
    RHO (pol)	 Te (eV)       Ne (1/m3)  Vtor_I (rad/s)        Ti1 (eV)     Ni1 (1/m3)	     Ni2 (1/m3)	     Ni3 (1/m3)
    """
    def __init__(self,dir_name, nrho, nion, A, Z):
        profiles.__init__(self)
        self.flag_ni2 = 1
        self.flag_ni3 = 1
        self.A = A
        self.Z = Z
        self.nrho = nrho
        self.nion = nion
        
        #self.ni_in = np.zeros((self.nion, 2, self.nrho), dtype=float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)

        te_fname = dir_name+'/te.dat'
        self.te_in=np.loadtxt(te_fname, unpack=True)
        nrho_in = np.size(self.te_in[0,:])
        self.ni_in=np.zeros((self.nion, 2,nrho_in), dtype=float)

        ne_fname = dir_name+'/ne.dat'
        vt_fname = dir_name+'/vtor.dat'
        ti_fname = dir_name+'/ti.dat'        
        i_fname = dir_name+'/ni1.dat'
        for i in range(self.nion):
            i_fname = dir_name+'/ni'+str(i+1)+'.dat'
            tmp_ni = np.loadtxt(i_fname, unpack=True)
            if np.shape(tmp_ni)[1]<nrho_in:
                tmp_ni2 = np.zeros((2,nrho_in))
                tmp_ni2[1,0:np.shape(tmp_ni)[1]] = tmp_ni[1,:]
                tmp_ni2[0,:] =np.linspace(0,1.2, nrho_in)
                tmp_ni = tmp_ni2
            self.ni_in[i,:,:] = tmp_ni
            
        self.rho=np.linspace(0,1,num=self.nrho)
        self.te_in=np.loadtxt(te_fname, unpack=True)
        self.ne_in=np.loadtxt(ne_fname, unpack=True)
        self.vt_in=np.loadtxt(vt_fname, unpack=True)
        self.ti_in=np.loadtxt(ti_fname, unpack=True)
        #self.ni_in[0,:,:] = np.loadtxt(i_fname, unpack=True)
        self.coll_mode = np.ones(self.nion+1)
       
        self.smooth()
   
    def smooth(self):
        """
        smooth input data to data wanted
        """
    
        self.te=self._spline(self.te_in[0,:], self.te_in[1,:], self.rho)
        self.ne=self._spline(self.ne_in[0,:], self.ne_in[1,:], self.rho)
        self.ti=self._spline(self.ti_in[0,:], self.ti_in[1,:], self.rho)
        self.vt=self._spline(self.vt_in[0,:], self.vt_in[1,:], self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.ni_in[i,0,:], self.ni_in[i,1,:], self.rho)
        self._extrapolate()


class matlab_profiles(profiles):
    """
    Function to write the profile file for ASCOT from a matlab file
    Pietro reads metis output and produces the matlab files that should be read here.
    
    INPUT:

    
    # Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation
    # range must cover [0,1] of normalised poloidal rho. It can exceed 1. 
    # 18Jan08 for testing (first 3 lines are comment lines)
    111	3	# Nrad,Nion
    1 1 6		# ion Znum				
    1 2 12		# ion Amass
    1 1 1 1		# collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons
    RHO (pol)	 Te (eV)       Ne (1/m3)  Vtor_I (rad/s)        Ti1 (eV)     Ni1 (1/m3)	     Ni2 (1/m3)	     Ni3 (1/m3)
    """
    def __init__(self,inf_name, nrho, convert):
        
        profiles.__init__(self)
        self.nrho = nrho
        self.rho = np.linspace(0,1,nrho, dtype=float)
        infile = sio.loadmat(inf_name)
        plasma = infile['plasma']
        p1d = plasma['p1d'][0,0]
        self.Z = plasma['znum'][0,0][:,0]
        self.A = plasma['anum'][0,0][:,0]
        self.coll_mode = plasma['colls'][0,0][:,0]
        self.nion = len(self.Z)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)

        self.rho_in = p1d['rho'][0,0][:,0]
        self.ni_in = np.zeros((self.nion, len(self.rho_in)),dtype=float)

        self.te_in = p1d['te'][0,0][:,0]
        self.ne_in = p1d['ne'][0,0][:,0]
        self.vt_in = np.zeros(len(self.rho_in))
        print "VTOR SET TO 0!"
        self.ti_in = p1d['ti'][0,0][:,0]
        for i in range(self.nion):
            self.ni_in[i, :] = p1d['ni'][0,0][:,i]
        self.smooth()

    def smooth(self):
        """
        smooth input data to data wanted
        """
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)
        self._extrapolate()
