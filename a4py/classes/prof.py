"""
matteo.vallar@igi.cnr.it - 02/2019

Class for profiles I/O for ascot&bbnbi
"""
from __future__ import print_function
import numpy as np
import h5py, time, platform
from scipy import interpolate, integrate
import scipy.io as sio
import a4py.classes.Bfield as ascot_Bfield
import a4py.plot.input.plot_input as plot_input

# if on TCV, load also mdsplus
cluster = platform.uname()[1]
if cluster[-7:] == 'epfl.ch':
    import MDSplus as mds

colours = ['k','b','r','c','g',\
           'k','b','r','c','g',\
           'k','b','r','c','g']

class profiles:
    """Profiles superclass:

    This class contains the methods shared among the other different classes
    
    Parameters:
        None

    Attributes:
        |  rho (array):
        |  te (array):
        |  ti (array):
        |  ne (array):
        |  ni (array):
        |  ni1 (array):
        |  ni2 (array):
        |  ni3 (array):
        |  vtor (array):
        |  zeff (array):
        |  nion (array):
        |  Z (array):
        |  A (array):
        |  collmode (array):   

    Methods:
        |  write_input: write input.plasma_1d file
        |  plot_profiles: Plot the profiles

    """
    def __init__(self):
        """Initialisation of profiles superclass

        Initialises the quantities used in the class

        """
        self.rho=[]
        self.te=[]
        self.ti=[]
        self.ne=[]
        self.ni=[]
        self.ni1=[]
        self.ni2=[]
        self.ni3=[]
        self.vtor=[]
        self.zeff=[]

        self.nion=1
        self.Z=[]
        self.A=[]
        self.coll_mode=[]
      

    def write_input(self, suffix=''):
        """ write input.plasma_1d file

        Method to write the input.plasma_1d file.
        output file would be input.plasma_1d
        
        Parameters:
            None

        Returns:
            None
            
        """
        
        out_fname = "input.plasma_1d"+suffix
        with open(out_fname, 'w+') as outfile:
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
            print(data)
            print("if i don't print, it won't work")
            np.savetxt(outfile, data, fmt='%.5e')

    def compute_average(self):
        """compute averages
        
    
        Parameters:
            None

        Returns:
            None
        """
        ind = self.rho<=1.
        vol = self._spline(self._rho_vol, self._volumes, self.rho[ind])
        ne_avg = np.trapz(self.ne[ind]*vol)
        ne_avg /= np.sum(vol)
        te_avg = np.trapz(self.te[ind]*vol)
        te_avg /= np.sum(vol)

        print('average ne (10e20 m-3): ', ne_avg*1.e-20)
        print('average te (      keV): ', te_avg*1.e-3)
        

    def _plot_time(self, tit):
        """function to plot with title

        This function uses the plot_profiles function to plot the profiles with a title

        Parameters:
            tit (str): title wanted on the plot (usually time)

        Returns:
            None

        """
        tit=str(tit)
        self.plot_profiles()

    def plot_profiles(self, fig=0):
        """Plot the profiles
        
        This function makes a plot with ne, Te, Ti, ni(eventually nimp) on 4 different frames

        Parameters:
            |  f (object): the plt.figure object where to plot (useful to overplot). Undefined is initialized to 0 and it means to do a new figure
            |  title (str): title of the figure
        Return:
            None

        """
        plot_input.plot_profiles(self, fig)


    def _spline(self,rho,data, rho_new):
        """ Splines of input

        Private method to evaluate splines of the data input, in order to make an output in the right rho grid

        Parameters:
            |  rho (array): x array where the data are defined
            |  data (array): data to be spline-d
            |  rho_new (array): new x array where data are wanted
        Returns:
            data_new (array): new data array desired

        """
        dummy    = interpolate.InterpolatedUnivariateSpline(rho, data, ext=0)
        data_new = dummy(rho_new)
        return data_new

    def _extrapolate(self):
        """ extrapolate over rho=1

        Private method that does an exponential fit outside rho=1 surface.
        The decay length is set to 0.01: y_out = y(rho=1)*exp(-(x-1.)/dec_l).
        It does it for all the quantities (ne, ni, nimp, te, ti, vtor).
        It's better to give it as input instead of letting it to ASCOT.

        Parameters:
            None
        Returns:
            None

        """
        x = np.linspace(1.001, 1.2, self.nrho/5)
        rho1 = self.rho # rho up to 1
        dec_l = 0.01
        ni_ov = np.zeros((self.nion, len(x)), dtype=float)
        ninew = np.zeros((self.nion, self.nrho+len(x)),dtype=float)
        ne_ov1 = self.ne[self.nrho-1]*np.exp(-((x-1.)/dec_l))
        te_ov1 = self.te[self.nrho-1]*np.exp(-(x-1.)/dec_l)
        ti_ov1 = self.ti[self.nrho-1]*np.exp(-(x-1.)/dec_l)
        vt_ov1 = self.vt[self.nrho-1]*np.exp(-(x-1.)/dec_l)
        for i in range(self.nion):
            ni_ov[i,:] = self.ni[i,self.nrho-1]*np.exp(-(x-1.)/dec_l)
            ninew[i,:] = np.concatenate([self.ni[i,:], ni_ov[i,:]])
        self.ni = ninew
        self.rho = np.concatenate([rho1, x])
        self.nrho = len(rho1)+len(x)
        self.ne  = np.concatenate([self.ne, ne_ov1])
        self.te  = np.concatenate([self.te, te_ov1])
        self.ti  = np.concatenate([self.ti, ti_ov1])
        self.vt  = np.concatenate([self.vt, vt_ov1])

    def _ion_densities(self):
        """Computes C and D ion densities

        Compute C ion densities starting from ni_in, ne_in and zeff.
        |  Solving the following system (valid only if D and C(6+) are the ion species):
            (1) Zeff = sum(ni*Zi**2)/sum(ni*Zi)
            (2) ne   = nD + 6*nC
        
        Parameters:
            None
        Arguments:
            None

        """
        nD = self.ne_in*(6-self.zeff_in)/(5.)
        nC = self.ne_in*(self.zeff_in-1)/(30.)
        nC[np.where(nC<0)]=0.
        print("nC/nD: "+str(np.mean(nC/nD)*100.)+" %")
        self.ni_in[0,:] = nD
        self.ni_in[1,:] = nC

class h5_profiles(profiles):   
    """ Class for handling the profiles data from h5 file

    This class reads the profiles from a h5 file ascot-like.

    Bases:
        profiles
    Parameters:
        |  infile_name (str): name of the h5 file to use
        |  nrho (int): number of rho points to use as output (up to 1, outside rho is 1/5 of nrho)
    Attributes:
        |  Inherited from profiles
        |  ne_in (array): input ne
        |  te_in (array): input te
        |  ti_in (array): input ti
        |  ni_in (array): input ni
        |  vt_in (array): input vtor (if present, otherwise 0)
        |  zeff_in (array): input zeff (if present, otherwise 0)
        
    Methods:
        |  read_h5: reads h5 file 
        |  smooth: calls smoothing routines from profiles for h5 input data

    Notes:
        DATA in h5 file (09/01/2017, ascot4)

        |  /plasma                  Group
        |  /plasma/1d               Group
        |  /plasma/1d/ne            Dataset {1373}
        |  /plasma/1d/ni            Dataset {1373, 3}
        |  /plasma/1d/rho           Dataset {1373}
        |  /plasma/1d/te            Dataset {1373}
        |  /plasma/1d/ti            Dataset {1373}
        |  /plasma/1d/vtor          Dataset {1373}
        |  /plasma/1d/zeff          Dataset {1373}
        |  /plasma/anum             Dataset {3}
        |  /plasma/colls            Dataset {4}
        |  /plasma/znum             Dataset {3}

    """

    def __init__(self, infile_name, nrho):
        profiles.__init__(self)
        #defining dictionary for handling data from h5file
        self.labdict = {"rho":"/plasma/1d/rho",\
                        "ne":"/plasma/1d/ne","ni":"/plasma/1d/ni",\
                        "te":"/plasma/1d/te","ti":"/plasma/1d/ti",\
                        "vtor":"/plasma/1d/vtor","zeff":"/plasma/1d/zeff",\
                        "a_ions":"/plasma/anum", "znum":"/plasma/znum"\
                        }
        self.inf_name = infile_name
        self.nrho = nrho

        try:
            self.read_h5()
        except IOError:
            "Not h5 file given as input"

    def read_h5(self):
        """ Reads h5
        
        This method reads the profiles from ascot.h5 file
        
        Parameters:
            None
        Returns:
            None

        Note:
            It calls the spline function at the end

        """
        infile = h5py.File(self.inf_name,'r')

        vardict = self.labdict
        #store data with the correct labels
        for k in infile['plasma/1d'].keys():
            try:
                vardict[k] = infile[self.labdict[k]].value
            except:
                vardict[k] = []

        vardict['a_ions']=infile['/plasma/anum'].value
        vardict['znum']=infile['/plasma/znum'].value
        

        self.rho_in = vardict['rho']
        self._rho_vol = infile['distributions/rhoDist/abscissae/dim1'].value[1:]
        self._volumes = infile['distributions/rhoDist/shellVolume'].value
        self.nrho_in = np.size(self.rho_in)

        if vardict['a_ions'][0]!='/':
            self.nspec = len(vardict['a_ions'])
        else:
            self.nspec = vardict['ni'].shape[1]
        print("Number of ions: ", self.nspec)
        if len(vardict['a_ions'])!=len(vardict['znum']):
            print("ERROR! array of A and Z don't have the same length")

        self.A = vardict['a_ions']
        self.Z = vardict['znum']
        self.nion = self.nspec
        
        self.te_in  = vardict['te'][:]
        self.ne_in  = vardict['ne'][:]
        self.ti_in  = vardict['ti'][:]
        ni1_in  = vardict['ni'][:,0]
        self.ni_in = np.zeros((self.nion, self.nrho_in),dtype=float)
        self.ni_in[0,:] = ni1_in
        if self.nion==2:
            ni2_in  = vardict['ni'][:,1]
            self.ni_in[1,:] = ni2_in
        elif self.nion==3:
            ni2_in  = vardict['ni'][:,1]
            ni3_in  = vardict['ni'][:,2]
            self.ni_in[1,:] = ni2_in
            self.ni_in[2,:] = ni3_in

        try:
            self.vt_in  = vardict['vtor']
        except:
            self.vt_in    = np.zeros(self.nrho_in,dtype=float)

        try:
            self.zeff_in  = vardict['zeff'][:]
        except:
            self.zeff_in  = np.zeros(self.nrho_in,dtype=float)

        self.ni = np.zeros((self.nion, self.nrho),dtype = float)
        self.spline()

    def spline(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        self.rho = np.linspace(0,1,self.nrho)
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)
        try:
            self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        except:
            self.vt = np.zeros(self.nrho, dtype=float)
        self.zeff = self._spline(self.rho_in, self.zeff_in, self.rho)

        self._extrapolate()


class dat_profiles(profiles):
    """ reads dat profile
    
    Class to handle the profiles for ASCOT from a series of
    ascii file in the format (rho, quantity)
    
    Base:
        profiles

    Parameters:
        | dir_name (str): directory where to find the files
        |  nrho (int): number of rho points to use as output (up to 1, outside rho is 1/5 of nrho)
        |  nion (int): number of ions to use
        |  A (array): atomic masses of ions
        |  Z (array): charge of ions
    
    Arguments:
        None

    Note:
        |  The files in the folder should be stored as:
        |  te.dat, ne.dat, ti.dat, ni1.dat, ni2.dat,...

    """
    def __init__(self,dir_name, nrho, nion, A, Z):
        profiles.__init__(self)
        self.flag_ni2 = 1
        self.flag_ni3 = 1
        self.A = A
        self.Z = Z
        self.nrho = nrho
        self.nion = nion
        
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
            print(np.shape(tmp_ni), nrho_in, tmp_ni)
            if np.shape(tmp_ni)[1]<nrho_in:
                 tmp_ni2 = np.zeros((2,nrho_in))
                 tmp_ni2[1,0:np.shape(tmp_ni)[1]] = tmp_ni[1,:]
                 tmp_ni2[0,:] =np.linspace(0,1.2, nrho_in)
                 tmp_ni = tmp_ni2
            self.ni_in[i,:,:] = tmp_ni
            
        self.rho=np.linspace(0,1,num=self.nrho)
        self.ne_in=np.loadtxt(ne_fname, unpack=True)
        self.vt_in=np.loadtxt(vt_fname, unpack=True)
        self.ti_in=np.loadtxt(ti_fname, unpack=True)
        self.coll_mode = np.ones(self.nion+1)
       
        self.smooth()
   
    def smooth(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        self.te=self._spline(self.te_in[0,:], self.te_in[1,:], self.rho)
        self.ne=self._spline(self.ne_in[0,:], self.ne_in[1,:], self.rho)
        self.ti=self._spline(self.ti_in[0,:], self.ti_in[1,:], self.rho)
        self.vt=self._spline(self.vt_in[0,:], self.vt_in[1,:], self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.ni_in[i,0,:], self.ni_in[i,1,:], self.rho)
        self._extrapolate()


class matlab_profiles(profiles):
    """ handles matlab profiles from Pietro

    Function to write the profile file for ASCOT from a matlab file
    Pietro (pietro.vincenzi@igi.cnr.it) reads metis output and produces the matlab files that should be read here.

    Base:
        profiles

    Parameters:
        |  inf_name (str): directory where to find the files
        |  nrho (int): number of rho points to use as output (up to 1, outside rho is 1/5 of nrho)
    Arguments:
        None
  
    """

    def __init__(self,inf_name, nrho):
        
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
        print("VTOR SET TO 0!")
        self.ti_in = p1d['ti'][0,0][:,0]
        for i in range(self.nion):
            self.ni_in[i, :] = p1d['ni'][0,0][:,i]
        self.smooth()

    def smooth(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)
        self._extrapolate()


class ufiles(profiles):
    """
    Function to write the profile file for ASCOT from the Ufiles

    __init__(self,pre,sub,shot,nrho): correctly initalise and get the data from the Ufile s
    smooth(self):smooth input data to grid wanted
    """

    def __init__(self, shot, nrho, eqdsk_fname='', **kwargs):
        profiles.__init__(self)
        if eqdsk_fname=='':
            print("Needs an eqdsk to do conversion rhotor-rhopol")
            raise KeyError
        else:
            self.eqdsk_fname = eqdsk_fname
            self._phi2psi()
        self.nrho=nrho; self.shot=shot
        print("Opening values for shot ", self.shot)
        print("Always considering Carbon as impurity")
        #GETTING number of rho
        fname = 'N'+str(shot)+'.ELE'
        tmp_uf = uf.RU(fname)
        self.nrho_in = tmp_uf.nf
        self.rho_in = tmp_uf.values['Y']
        t_in   = tmp_uf.values['X']

        self.nion = 2 # D and C always included
        self.rho   = np.linspace(0,1,self.nrho, dtype = float)
        self.ne_in = np.zeros(self.nrho_in, dtype = float)
        self.ni_in = np.zeros((self.nion,self.nrho_in), dtype = float)
        self.te_in = np.zeros(self.nrho_in, dtype = float)
        self.ti_in = np.zeros(self.nrho_in, dtype = float)
        self.vt_in = np.zeros(self.nrho_in, dtype = float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        if 'ZEFF' in kwargs:
            self.zeff_in = np.full((self.nrho_in), kwargs['ZEFF'], dtype = float)
            print("Set Zeff from argument ZEFF: "+str(kwargs['ZEFF']))
        else:
            self.zeff_in = np.zeros(self.nrho_in, dtype = float)
           
        self.A = [2, 12]
        self.Z = [1, 6]
        self.coll_mode = np.ones(self.nion+1)
        pre_arr  = ['N', 'T', 'T', 'Z']
        suff_arr = ['ELE', 'ELE', 'ION', 'EFF']
        if 't' in kwargs:
            ind = np.argmin(t_in-kwargs['t'] < 0.)
            self.time = t_in[ind]
        else:
            if len(t_in)>1:
                print("More timeslices than 1, put kwargs with 't=time' to choose timeslice")
                raise ValueError
            else:
                ind = 0
                self.time = t_in[ind]
        for pre, suff in zip(pre_arr, suff_arr):
            fname = pre+str(shot)+'.'+str(suff)
            try:
                tmp_uf = uf.RU(fname)
            except:
                print("No File "+fname)
                continue
            if pre == 'N':
                self.ne_in = tmp_uf.fvalues[ind]*1e6 #conversion from cm^-3 to m^-3
            elif pre == 'T':
                if suff == 'ELE':
                    self.te_in = tmp_uf.fvalues[ind]
                else: #T ION
                    print(tmp_uf.fvalues)
                    if tmp_uf.fvalues.shape[1] == self.nrho_in:
                        self.ti_in = tmp_uf.fvalues[ind]
                    else:
                        self.ti_in = np.swapaxes(tmp_uf.fvalues, 0,1)[ind]
            elif pre == 'Z':
                self.zeff_in = tmp_uf.fvalues[0]
                print("Set Zeff from file "+fname)

        if np.mean(self.zeff_in)==0:
            self.zeff_in = np.full((self.nrho_in), 2., dtype = float)
            print("Set Zeff to 2")

        self.rhopol = self.param_psi(np.linspace(0,1,len(self.rho_in)))
        self.rhotor = self.rho_in
        self.rho_in = self.rhopol
       
        self._ion_densities()
        self._smooth()
    
    def _readeqdsk(self):
        """
        reads q and the poloidal flux from an eqdsk to convert phi2psi
        """
        dir_TCV = '/home/vallar/TCV/eqdsk/'
        try:
            b = ascot_Bfield.Bfield_eqdsk(dir_TCV+self.eqdsk_fname,129,129, 'TCV', COCOS=17)
            print("Opened ", dir_TCV+self.eqdsk_fname)
        except:
            print("Impossible to open ", self.eqdsk_fname)
            raise ValueError
        qprof_t   = np.abs(b.eqdsk.q)
        rho_eqdsk = b.eqdsk.rhopsi
        self.param_q = interpolate.interp1d(rho_eqdsk, qprof_t)

    def _phi2psi(self):
        """
        Converts psi 2 phi
        """
        try:
            self.param_q.mean()
        except:
            self._readeqdsk()
        tmpnum=100000
        locq   = self.param_q(np.linspace(0,1,tmpnum)) #augmenting precision near the core
        locphi = np.linspace(0,1,tmpnum)
        psi = integrate.cumtrapz(1/locq,locphi)
        psi = np.concatenate([[0], psi])
        psi = psi/max(psi)
        rhopsi = psi
        self.param_psi = interpolate.interp1d(np.linspace(0,1,tmpnum), rhopsi)

    def _smooth(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)
        self._extrapolate()


    def plot_time(self):
        title = 'Shot '+str(self.shot)+' | t = '+str(self.time)
        self._plot_time(title)

class SA_datfiles(profiles):
    """ class to handle SA datfiles, as read in repository

    This class is to read the files from CRONOS (usually produced by Jeronimo Garcia (jeronimo.garcia@cea.fr))
    
    Base:
        profiles
    
    Parameters:
        |  infile_name (str): name of the h5 file to use
        |  nrho (int): number of rho points to use as output (up to 1, outside rho is 1/5 of nrho)
        |  nion (int): number of ions to use
        |  A (array): atomic masses of ions
        |  Z (array): charge of ions

    Arguments:
        None

    Note
        rho(TOR)	ne	te	ti	zeff	psupra	nsupra	Jtot	jboot	jnbi	jec

    """
    def __init__(self, infile, nrho, nion, A, Z, shot, zeff=[0.]):
        profiles.__init__(self)

        self.A = A
        self.Z = Z
        self.nrho = nrho
        self.nion = nion
        #self.ni_in = np.zeros((self.nion, 2, self.nrho), dtype=float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        self.coll_mode = np.ones(self.nion, dtype=float)
        self.shot=shot
        #Read eqdsk to convert from rhotor to rhopol
        #self._readeqdsk(shot)

        lines = np.loadtxt(infile, skiprows=1, unpack=True)
        self.rho_in = lines[0,:]
        self.rhotor = self.rho_in
        #self.param_q = interpolate.interp1d(self.rhotor, lines[5,:])
        self._phi2psi()
        self.psipol = self.param_psi(self.rhotor) #This is psipol as function of rhotor
        self.rho_in = self.psipol**0.5 #this is rhopol as function of rhotor
        _rhop = interpolate.interp1d(self.rhotor, self.rho_in) #parameters for re-gridding
        self.rho = _rhop(np.linspace(0,1,self.nrho))
        
        self.ne_in  = lines[1,:]
        self.te_in  = lines[2,:]
        self.ti_in  = lines[3,:]
        if zeff[0]==0:
            self.zeff_in  = lines[4,:]
        else:
            if len(zeff)==1:
                self.zeff_in = np.full(len(self.rho_in), zeff, dtype=float)
            else:
                p=interpolate.interp1d(np.linspace(0,1,len(zeff)),zeff)
                self.zeff_in = p(self.rho_in)
        self.ni_in  = np.zeros((self.nion, len(self.rho_in)),dtype=float)
        self.vt_in = np.zeros(len(self.rho_in),dtype=float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        if len(self.Z)>1:
            self._ion_densities()
        self.smooth()

    def _readeqdsk(self, shot):
        """ Reads eqdsk

        reads q and the poloidal flux from an eqdsk to convert phi2psi

        Parameters:
            shot: identifier of scenario (for SA 002, 003, etc)
        Arguments:
            None
        """
        if shot==3 or shot==2:
            dir_JT = '/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/'
            eqdsk_fname = 'Equil_JT60_prova01_e_refined.eqdsk'
        elif shot==4:
            dir_JT = '/home/vallar/JT60-SA/004_2/input_from_EUrepository/'
            eqdsk_fname = 'JT-60SA_scenario4_uptowall.geq'
        elif shot==5:
            dir_JT = '/home/vallar/JT60-SA/005/input_from_EUrepository/'
            eqdsk_fname = 'JT-60SA_scenario5_eqdsk'
        try:
            b = ascot_Bfield.Bfield_eqdsk(dir_JT+eqdsk_fname,129,129, 'JT60SA', COCOS=3)
            print("Opened ", dir_JT+eqdsk_fname)
        except:
            print("Impossible to open ", eqdsk_fname)
            raise ValueError

        qprof_t   = b.eqdsk.q
        rho_eqdsk = b.eqdsk.rhopsi
        self.param_q = interpolate.interp1d(rho_eqdsk, qprof_t)

    def smooth(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)

        self.zeff = self._spline(self.rho_in, self.zeff_in, self.rho)

        self._extrapolate()

    def _phi2psi(self):
        """Converts psi 2 phi
        
        Converts in between coordinates using the following relation (phi=toroidal flux, psi=poloidal flux)
        psi = int(1/q, phi)

        Paramters:
            None
        Arguments:
            None
        """
        try:
            locq   = self.param_q(self.rhotor)
        except:
            self._readeqdsk(self.shot)
            locq   = self.param_q(self.rhotor)
            
        locphi = self.rhotor**2
        psi = integrate.cumtrapz(1/locq,locphi)
        psi = np.concatenate([[0], psi])
        psi = psi/max(psi)
        self.param_psi = interpolate.interp1d(self.rhotor, psi)   
        

        # tmpnum=100000
        # locq   = self.param_q(np.linspace(0,1,tmpnum)) #augmenting precision near the core
        # locphi = self.rhotor**2
        # locphi_p = interpolate.interp1d(np.linspace(0,1,len(locphi)),locphi)
        # locphi = locphi_p(np.linspace(0,1,tmpnum))
        # psi = integrate.cumtrapz(1/locq,locphi)
        # psi = np.concatenate([[0], psi])
        # psi = psi/max(psi)
        # rhopsi = psi**0.5
        # self.param_psi = interpolate.interp1d(np.linspace(0,1,tmpnum), rhopsi)

class SA_datfiles_datascenario(profiles):
    """ class to handle SA datfiles, as read in repository

    This class is to read the files from CRONOS (usually produced by Jeronimo Garcia (jeronimo.garcia@cea.fr))
    
    Base:
        profiles
    
    Parameters:
        |  infile_name (str): name of the h5 file to use
        |  nrho (int): number of rho points to use as output (up to 1, outside rho is 1/5 of nrho)
        |  nion (int): number of ions to use
        |  A (array): atomic masses of ions
        |  Z (array): charge of ions

    Arguments:
        None

    Note
    rho(TOR)	psi	ne	ni	te	ti	Jpar	q	"Pnbi,el"	"Pnbi,ion"	"Pnbi,supra"	"Nbi,nfast"	"J,NBI"	"Pecrh,el"	Jeccd
    ni is assumed to be the D density

    """
    def __init__(self, infile, nrho, shot):
        profiles.__init__(self)
        A=[2,12]; Z=[1,6]
        nion=2
        print("Using A=2,12 and Z=1,6")
        self.A = A
        self.Z = Z
        self.nrho = nrho
        self.nion = nion
        self.rho = np.linspace(0,1,self.nrho, dtype=float)
        #self.ni_in = np.zeros((self.nion, 2, self.nrho), dtype=float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        self.coll_mode = np.ones(self.nion, dtype=float)
        
        lines = np.loadtxt(infile, skiprows=1, unpack=True)
        self.rho_in = lines[0,:]
        self.rhotor = self.rho_in
        self.param_q = interpolate.interp1d(self.rhotor, lines[7,:])
        self._phi2psi()
        self.psipol = self.param_psi(np.linspace(0,1,len(self.rhotor)))
        self.rho_in = self.psipol**0.5 #this is rhopol as function of rhotor
        _rhop = interpolate.interp1d(self.rhotor, self.rho_in) #parameters for re-gridding
        self.rho = _rhop(np.linspace(0,1,self.nrho))

        self.ne_in  = lines[2,:]
        self.te_in  = lines[4,:]
        self.ti_in  = lines[5,:]
        self.ni_in  = np.zeros((self.nion, len(self.rho_in)),dtype=float)
        self.ni_in[0,:] = lines[3,:]
        self.vt_in = np.zeros(len(self.rho_in),dtype=float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        if len(self.Z)>1:
            self._ion_densities_datafiles()
        self.smooth()


    def _ion_densities_datafiles(self):
        """Computes C and D ion densities
        Given the ne, nD, computes NC
        nC = (ne-nD)/6.
       
        Parameters:
            None
        Arguments:
            None

        """
        ne = self.ne_in
        nD = self.ni_in[0,:]
        nC = (ne-nD)/6.
        print("nC/nD: "+str(np.mean(nC/nD)*100.)+" %")
        self.ni_in[0,:] = nD
        self.ni_in[1,:] = nC

    def _readeqdsk(self, shot):
        """ Reads eqdsk

        reads q and the poloidal flux from an eqdsk to convert phi2psi

        Parameters:
            shot: identifier of scenario (for SA 002, 003, etc)
        Arguments:
            None
        """
        if shot==3 or shot==2:
            dir_JT = '/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/'
            eqdsk_fname = 'Equil_JT60_prova01_e_refined.eqdsk'
        elif shot==4:
            dir_JT = '/home/vallar/JT60-SA/004_2/input_from_EUrepository/'
            eqdsk_fname = 'JT-60SA_scenario4_uptowall.geq'
        elif shot==5:
            dir_JT = '/home/vallar/JT60-SA/005/input_from_EUrepository/'
            eqdsk_fname = 'JT-60SA_scenario5_eqdsk'
        try:
            b = ascot_Bfield.Bfield_eqdsk(dir_JT+eqdsk_fname,129,129, 'JT60SA', COCOS=3)
            print("Opened ", dir_JT+eqdsk_fname)
        except:
            print("Impossible to open ", eqdsk_fname)
            raise ValueError

        qprof_t   = b.eqdsk.q
        rho_eqdsk = b.eqdsk.rhopsi
        self.param_q = interpolate.interp1d(rho_eqdsk, qprof_t)

    def smooth(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        self.vt = self._spline(self.rho_in, self.vt_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)

        #self.zeff = self._spline(self.rho_in, self.zeff_in, self.rho)

        self._extrapolate()


    def _phi2psi(self):
        """Converts psi 2 phi
        
        Converts in between coordinates using the following relation (phi=toroidal flux, psi=poloidal flux)
        psi = int(1/q, phi)

        Paramters:
            None
        Arguments:
            None
        """
        try:
            locq   = self.param_q(self.rhotor)
        except:
            self._readeqdsk(self.shot)
            locq   = self.param_q(self.rhotor)
            
        locphi = self.rhotor**2
        psi = integrate.cumtrapz(1/locq,locphi)
        psi = np.concatenate([[0], psi])
        psi = psi/max(psi)
        self.param_psi = interpolate.interp1d(self.rhotor, psi)   

        #tmpnum=100000
        # locq   = self.param_q(np.linspace(0,1,tmpnum)) #augmenting precision near the core
        # locphi = np.linspace(0,1,tmpnum)
        # psi = integrate.cumtrapz(1/locq,locphi)
        # psi = np.concatenate([[0], psi])
        # psi = psi/max(psi)
        # rhopsi = psi**0.5
        # self.param_psi = interpolate.interp1d(np.linspace(0,1,tmpnum), rhopsi)


class TCV_datfiles(profiles):
    """
    Reads a dat file with the columns as follows:
    rho pol, ne [m^-3], te[eV], ti[eV]
    """
    def __init__(self, infile, **kwargs):
        profiles.__init__(self)
        
        self.A = [2, 12]
        self.Z = [1, 6]
        self.nion = 2
        #self.ni_in = np.zeros((self.nion, 2, self.nrho), dtype=float)

        
        lines = np.loadtxt(infile, skiprows=1, unpack=True)
        self.rho_in = lines[0,:] #already poloidal rho
        self.nrho_in = len(self.rho_in)
        self.rho = self.rho_in
        self.nrho = self.nrho_in
        self.ne_in  = lines[1,:]
        self.te_in  = lines[2,:]
        self.ti_in  = lines[3,:]

        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        self.vt = np.zeros((self.nrho), dtype=float)
        self.coll_mode = np.ones(self.nion, dtype=float)
        
        if 'ZEFF' in kwargs:
            self.zeff_in = np.full((self.nrho_in), kwargs['ZEFF'], dtype = float)
            print("Set Zeff from argument ZEFF: "+str(kwargs['ZEFF']))
        else:
            print('Zeff not set, putting it to 0')
            self.zeff_in = np.zeros(self.nrho_in, dtype = float)

        self.ni_in  = np.zeros((self.nion, len(self.rho_in)),dtype=float)
        self.ni = np.zeros((self.nion, self.nrho), dtype=float)
        if len(self.Z)>1:
            self._ion_densities()
        self.smooth()

    def smooth(self):
        """ spline input data to grid wanted
        
        For each variable the input array is splined and put in output to
        desired grid. This is specific for h5 files

        Parameters:
            None
        Attributes:
            None
        
        Note:
            The _extrapolate private method is called
        
        """
        self.te = self._spline(self.rho_in, self.te_in, self.rho)
        self.ne = self._spline(self.rho_in, self.ne_in, self.rho)
        self.ti = self._spline(self.rho_in, self.ti_in, self.rho)
        for i in range(self.nion):
            self.ni[i,:]=self._spline(self.rho_in, self.ni_in[i,:], self.rho)

        self.zeff = self._spline(self.rho_in, self.zeff_in, self.rho)

        self._extrapolate()


class TCV_mds(profiles):
    """ Reads profiles from TCV tree

    Connects to TCV tree and reads the data from there. ne, Te from thomson, Ti from cxrs

    Parameters:
        indict (dict): 
            |  'shot' (int): shot number,
            |  't' (float): time
            |  'nrho' (int): number of rho points to use as output (up to 1, outside rho is 1/5 of nrho)
            |  'zeff' (float, optional): zeff in that time instant
    Arguments:
        None
    Notes:
        |  'ne':'tcv_shot::top.results.thomson.profiles.auto:ne'
        |  'te':'tcv_shot::top.results.thomson.profiles.auto:te'
        |  'ti':'tcv_shot::top.results.cxrs.proffit:ti'
    
    """
    def __init__(self, indict):
        self.indict = indict
        self.shot = indict['shot']
        self.t = indict['t']
        self.nrho = indict['nrho']
        self.rho = np.linspace(0, 1, num=self.nrho)
        try:
            self.zeff = indict['zeff']
        except:
            print("No Zeff set! setting it to 2.")
            self.zeff = 2.
        self.nion = 2
        self.A = [2, 12 ]
        self.Z = [1, 6  ]
        self.coll_mode = np.ones(self.nion+1)
        # open the tree
        self.tree = mds.Tree('tcv_shot', self.shot)

        # we build the appropriate dictionary similarly to what done
        # for the 1D signal
        self.signals = {'ne': {'string': r'\tcv_shot::top.results.thomson.profiles.auto:ne'},
                        'te': {'string': r'\tcv_shot::top.results.thomson.profiles.auto:te'},
                        'ti': {'string': r'\tcv_shot::top.results.cxrs.proffit:ti'}}

        print("\n")
        print("===================")
        print("Initialisation of 2D signals  Done")
        print("===================")
        print("\n")


    def _getBivecSpline(self):
        """ reads and spline of input signals
        Private method for reading the signal storing their bivacspline
        representation on a grid (time, rho_tor)

        Parameters:
            None
        Arguments:
            None
        """
        self._brep = {}

        for k in self.signals.keys():
            print('Reading signal ' + self.signals[k]['string'])
            tim = self.tree.getNode(self.signals[k]['string']).getDimensionAt(0).data()
            if k=='ti':
                tim = self.tree.getNode(self.signals[k]['string']).getDimensionAt(1).data()
            _idx = np.argmin(tim-self.t < 0)
            tim = tim[_idx]
            if k=='ti':
                # CXRS
                data = self.tree.getNode(self.signals[k]['string']).data()[_idx,:]
                rhop = self.tree.getNode(self.signals[k]['string']).getDimensionAt(0).data()[_idx,:]
                indsort = np.argsort(rhop)
                data = data[indsort]
                rhop = rhop[indsort]
            else:
                data = self.tree.getNode(self.signals[k]['string']).data()[:, _idx]
                rhop = self.tree.getNode(r'\results::thomson.profiles.auto:rho').data()
            dummy = interpolate.interp1d(rhop, data, fill_value='extrapolate')
            self._brep[k] = dict([('spline', dummy)])

            
    def read_2d(self):
        """ read 2d signals

        Method to get the signal defined in an attribute of the
        class with the appropriate resolution in rho_toroidal
        and time. It create the attribute self.rsig with a dictionary
        with all the signals with corresponding time and rho basis

        Parameters:
            None
        Arguments:
            None
        """
        try:
            self._brep
        except:
            self._getBivecSpline()
            
        self.rsig = {}
        for k in self.signals.keys():
            y = self._brep[k]['spline'](self.rho)
            self.rsig[k] = dict([('signal',y), 
                                 ('rho', self.rho)])

        self._tosuperclass()
            
        print("\n")
        print("===================")
        print("END READING 2D")
        print("===================")
        print("\n")


    def _tosuperclass(self):
        """interface between MDS and superclass

        Converts the input read from MDS to superclass useful input

        Parameters:
            None
        Arguments:
            None
        """        
        self.ne_in = self.rsig['ne']['signal']
        self.ne = self.ne_in
        self.te_in = self.rsig['te']['signal']
        self.ti_in = self.rsig['ti']['signal']
        self.ni_in = np.zeros((self.nion, len(self.ne_in)),dtype=float)
        self.zeff_in = np.full(self.nrho, self.zeff)
        self.vt_in = np.zeros(len(self.ne_in),dtype=float)
        self.vt = np.zeros(len(self.ne_in),dtype=float)
        self._ion_densities()
        self.ni = self.ni_in
        self.te = self.te_in
        self.ti = self.ti_in
        
        # no need to smooth since they are already smoothed
        self._extrapolate()
                        
    def write_tcv(self):
        """ write output for tcv
        
        Uses the write_input routine to write output with the time used
        
        Parameters:
            None
        Arguments:
            None
        """
        suffix = '_'+str(self.shot)+'_'+str(int(self.t*1e3))
        self.write_input(suffix=suffix)

class input_datfiles(profiles):
    """
    Reads an input.plasma_1d file:
    rho pol, ne [m^-3], te[eV], ti[eV], ni [m^-3]
    """
    def __init__(self,infile):
        """
        """
        self.read_input(infile)
        
    def read_input(self,infile):    
        fin = open(infile,"r")
        lines=fin.readlines()

        #First three lines are commented
        # fourth line has {:d}\t{:1d}\t# Nrad,Nion \n'.format(self.nrho,self.nion)
        self.nrho = int(lines[3].split()[0])
        self.nion = int(lines[3].split()[1])

        # fifth line has '\t\t# ion Znum \n'
        # sixth line has '\t\t# ion Amass \n'
        # seventh line has '# collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons \n'

        self.Z = np.zeros(self.nion)
        self.A = np.zeros(self.nion)
        for i in range(self.nion):
            self.Z[i] = float(lines[4].split()[i])
            self.A[i] = float(lines[5].split()[i])

        #eight line has RHO (pol)      Te (eV)        Ne (1/m3)      Vtor_I (rad/s) Ti1 (eV)       Ni1 (1/m3)     Ni2 (1/m3)
        data=np.zeros((self.nrho, 5+self.nion), dtype=float)
        for i,el in enumerate(lines[8:]):
            tmp = el.split()
            data[i,:] = [float(tmp[j]) for j in range(len(tmp))]
        data=data.T
        self.rho = data[0,:]
        self.te = data[1,:]; self.ne = data[2,:]
        self.vtor = data[3,:]
        self.ti = data[4,:]
        self.ni = np.zeros((self.nrho, self.nion))
        self.ni[:,0] = data[5,:]
        for i in range(self.nion-1):
            self.ni[:,i+1] = data[i+6,:]
        self.ni = np.transpose(self.ni)
