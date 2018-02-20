"""
matteo.vallar@igi.cnr.it - 11/2017

Class for distributions
two classes inside:
distribution_1d(h5 ascot file)
distribution_2d(h5 ascot file with 2D distributions in it)

"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import ticker, colors
from scipy import interpolate
import os.path, math, time
import collections

import ascot_particles

colours = ['k','g','c','b','r']

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
                  
class distribution_1d:
    """
    Class for handling the <distributions> data

    METHODS:
    __init__(self, infile): getting basic variables, checking if rhoDist 
        in hdf5 file
    rhodists(self): Method to get the data from the ascot file
    TCV_calc_FIBP(self, plot_flag, *args): Compute the FIBP in TCV case. 
        Plot_flag=1 to plot, if given as *args shot and run it plots 
        the FIBP from a different BBNBI h5 file
    TCV_plot_all(self): Plot quantities w/o different ion species (i.e. j,p,ecc.)


    plot_current(self): Plot the induced beam current density
    plot_power(self): Plot just the deposited power density
    plot_torque(self): Plot just the deposited torque density
    plot_totalcurrent(self): Plot sum of the induced beam current density
    plot_totalpower(self): Plot sum of the deposited power density
    plot_totaltorque(self): Plot sum of the torque density

    store_groupdis_to_ascii(self, fname): Store total data distribution to ASCII

    print_scalars(self): Print the scalars

    calc_eff(self): calculates the efficiency in the current drive
    
    HIDDEN METHODS:
    _evaluate_shellVol(self, rho_new): Evaluate the volumes at rho_new
    _evaluate_shellArea(self, rho_new): Evaluate the areas at rho_new 
    _group_beams(self): Group data distribution forgetting about beam index
    _calculate_scalar(self): Calculate scalar quantities



    DATA in h5 file (09/01/2017, ascot4)
    state can be both inistate or endstate

    /distributions/rhoDist   Group
    /distributions/rhoDist/abscissae Group
    /distributions/rhoDist/abscissae/dim1 Dataset {201}
    /distributions/rhoDist/abscissae/dim2 Dataset {2}
    /distributions/rhoDist/abscissae/dim3 Dataset {2}
    
    # name and unit for the /distributions/rhoDist/ordinates/
    # got with the command infile['/distributions/rhoDist/ordinates/'].keys()
    # for i in range():
    #     print infile['/distributions/rhoDist/ordinates/name_0000'+str(i)].value
    1  density m^{-3}
    2  energy density J m^{-3}
    3  parallel current density A/m^2
    4  toroidal current density A/m^2
    5  jxB Torque N m^{-2}
    6  jxB Torque, from ini/endstate N m^{-2}
    7  CX ion source m^{-3}
    8  CX ion energy source J m^{-3}
    9  CX neutral source m^{-3}
    10 CX neutral energy source J m^{-3}
    11 Pphi torque N m^{-2}
    12 parallel energy density J m^{-3}
    13 total toroidal current density A/m^2
    14 total pressure N m^{-2}
    15 parallel pressure N m^{-2}
    16 perpendicular pressure N m^{-2}
    17 finite Larmor radius torque N m^{-2}
    18 Thermalized particle density m^{-3}
    19 Thermalized particle energy density J m^{-3}
    20 Thermalized particle torque N m^{-2}
    21 Absorbed ICRH power W m^{-3}
    22 J.B
    23 Power deposition to electrons W m^{-3}
    24 Power deposition to background species  1 W m^{-3}
    25 Power deposition to background species  2 W m^{-3}
    26 Power deposition to background species  3 W m^{-3}
    27 Collisional torque deposition to electrons N m^{-2}
    28 Collisional torque deposition to background species  1 N m^{-2}
    29 Collisional torque deposition to background species  2 N m^{-2}
    30 Collisional torque deposition to background species  3 N m^{-2}


    # the groups here below are for debugging
    
    /distributions/accdist   Group
    /distributions/accdist/abscissae Group
    /distributions/accdist/abscissae/dim1 Dataset {101}
    /distributions/accdist/abscissae/dim2 Dataset {2}
    /distributions/accdist/abscissae/dim3 Dataset {2}
    /distributions/accdist/ordinate Dataset {1, 1, 1, 1, 1, 100, 1}
    /distributions/accdist/ordinates Group
    /distributions/accdist/ordinates/name_000001 Dataset {SCALAR}
    /distributions/accdist/ordinates/unit_000001 Dataset {SCALAR}
    
    /distributions/rhodtcolldist Group
    /distributions/rhodtcolldist/abscissae Group
    /distributions/rhodtcolldist/abscissae/dim1 Dataset {201}
    /distributions/rhodtcolldist/abscissae/dim2 Dataset {51}
    /distributions/rhodtcolldist/abscissae/dim3 Dataset {2}
    /distributions/rhodtcolldist/abscissae/dim4 Dataset {2}
    /distributions/rhodtcolldist/ordinate Dataset {1, 1, 1, 1, 50, 200, 1}
    /distributions/rhodtcolldist/ordinates Group
    /distributions/rhodtcolldist/ordinates/name_000001 Dataset {SCALAR}
    /distributions/rhodtcolldist/ordinates/unit_000001 Dataset {SCALAR}
    
    /distributions/rhodtdist Group
    /distributions/rhodtdist/abscissae Group
    /distributions/rhodtdist/abscissae/dim1 Dataset {201}
    /distributions/rhodtdist/abscissae/dim2 Dataset {51}
    /distributions/rhodtdist/abscissae/dim3 Dataset {2}
    /distributions/rhodtdist/abscissae/dim4 Dataset {2}
    /distributions/rhodtdist/ordinate Dataset {1, 1, 1, 1, 50, 200, 1}
    /distributions/rhodtdist/ordinates Group
    /distributions/rhodtdist/ordinates/name_000001 Dataset {SCALAR}
    /distributions/rhodtdist/ordinates/unit_000001 Dataset {SCALAR}
     """

    def __init__(self, infile_n):
        """
        getting basic variables, checking if rhoDist in hdf5 file
        """

        self.infile=h5py.File(infile_n)
        self.infile_n = infile_n
        rhonew=self.infile['plasma/1d/rho'][:]
        rhonew = rhonew[rhonew<1]
        try:
            self._volumes = self.infile['distributions/rhoDist/shellVolume'].value
            self._areas   = self.infile['distributions/rhoDist/shellArea'].value
            self._evaluate_shellVol(rhonew)
            self._evaluate_shellArea(rhonew)
        except:
            print("No /distributions/rhoDist/ in ", infile_n)    
        
        self.rhodists()
        self.fibp_particles = 0



    def rhodists(self):
        """
        Method to get the data from the ascot file
        """
        tree_path = '/distributions/rhoDist/'
        
        #This dictionary is usable if only one bulk species is present in the plasma
        self.name_dict = {'n':1, 'e_den':2,\
                          'jpar':3, 'jperp':4, \
                          'jxB':5, 'jxBstate':6, \
                          'CX_ionsource':7, 'CX_ionensource':8,\
                          'CX_neutsource':9, 'CX_neutensource':10,\
                          'torque':11,'par_e_den':12,\
                          'tot_tor_j':13, 'ptot':14, 'ppar':15, 'pperp':16,\
                          'flr_torque':17,\
                          'th_n':18, 'th_e_n':19, 'th_torque':20, 'abs_ICRH':21\
                          }
        self._check_dims()
        self.abscissae = {}
        self.abscissae = self.abscissae.fromkeys(self.infile['/distributions/rhoDist/abscissae'].keys(),0)
        for key in self.abscissae.keys():
            self.abscissae[key]=self.infile[tree_path+'/abscissae/'+str(key)].value
        self.rho = np.linspace(0,1,len(self.abscissae['dim1'])-1)
            
        ordinate = self.infile['/distributions/rhoDist/ordinate'].value
        
        #ADDING DIFFERENT ION SPECIES
        self.nions = len(self.infile['plasma/anum'][:])
        n_ions_more = self.nions-1
        self.name_dict['ctor_el']+= n_ions_more
        self.name_dict['ctor_i1']+= n_ions_more              
        for el in range(n_ions_more):
            k='pi'+str(el+2)
            self.name_dict[k]=24+el
            k2='ctor_i'+str(el+2)
            self.name_dict[k2]=26+el+1

        #self.slices structure WILL BECOME:
        #(injector, time, rho, type of distribution)
        self.slices = ordinate.reshape(ordinate.shape[-4], ordinate.shape[-3], ordinate.shape[-2], ordinate.shape[-1])
        self.n_inj = self.slices.shape[0]
        self.dim_num = len(self.infile['/distributions/rhoDist/ordinates/'].keys())
        self.dim_num = int(self.dim_num/2.) #this because names and units are doubled
        self.lab = np.array([], dtype='S32')
        self.uni = np.array([], dtype='S8')
        for i in range(self.dim_num):
            self.lab = np.append(self.lab, self.infile[tree_path+'ordinates/name_'+'{:06d}'.format(i+1)].value)
            self.uni = np.append(self.uni, self.infile[tree_path+'ordinates/unit_'+'{:06d}'.format(i+1)].value)
            
        #in infile['species/testParticle/origin'] there is an array with the ordered set of beams
        self._h5origins = self.infile['species/testParticle/origin'].value
        self._group_beams()


    def _check_dims(self):
        """
        In two different versions of ascot (e.g. 9215 and >9401) there is one
        rhodist more, which is J.B. Need to check for that to be the case
        """
        dimnum = len(self.infile['/distributions/rhoDist/ordinates/'].keys())
        dimnum -= len(self.infile['plasma/anum'][:])*2
        dimnum /= 2
        if dimnum==25:
            self.name_dict['pel'] = 22
            self.name_dict['pi1'] = 23
            self.name_dict['ctor_el'] = 24
            self.name_dict['ctor_i1'] = 25 # The last two lines are affected by the number of ion species
        else:
            self.name_dict['J.B'] = 22
            self.name_dict['pel'] = 23
            self.name_dict['pi1'] = 24
            self.name_dict['ctor_el'] = 25
            self.name_dict['ctor_i1'] = 26 # The last two lines are affected by the number of ion species
        return
       
       
    def plot_current(self):
        """
        Plot the induced beam current density
        """
        if self.n_inj==1:
            self.plot_totalcurrent()
        
    def plot_power(self):
        """
        Plot just the deposited power density
        """
        if self.n_inj==1:
            self.plot_totalpower() 
       
    def plot_torque(self):
        """
        Plot just the deposited torque density
        """
        if self.n_inj==1:
            self.plot_totaltorque()  

    def plot_totalcurrent(self):
        """
        Plot sum of the induced beam current density
        """
        i_tot = self.slices_summed[0,:,12]*1e-3
        if np.mean(i_tot)<0:
            i_tot = -1*i_tot
        plot_article(1,[self.rho, i_tot],[''],r'$\rho$', 'j (kA/$m^2$)', self.infile_n)

    def plot_totalpower(self):
        """
        Plot sum of the deposited power density
        """
        ind = self.name_dict['pel']-1
        if self.slices.shape[2] == 27:
            ind=ind-1
        pe = self.slices_summed[0,:,ind]*1e-3
        pi1 = self.slices_summed[0,:,ind+1]*1e-3
        if self.nions > 1:
            pi2 = self.slices_summed[0,:,ind+2]*1e-3
            if self.nions == 2:
                plot_article(3,[self.rho, pe, pi1, pi2],['el.', 'i1', 'i2'],r'$\rho$', 'p (kW/$m^3$)', self.infile_n)        
                
            elif self.nions == 3:
                pi3 = self.slices_summed[0,:,ind+3]*1e-3
                plot_article(4,[self.rho, pe, pi1, pi2, pi3],['el.', 'i1', 'i2', 'i3'],r'$\rho$', 'p (kW/$m^3$)', self.infile_n)        

        else:
            plot_article(2,[self.rho, pe, pi1],['el.', 'i1'], r'$\rho$', 'p (kW/$m^3$)', self.infile_n)        

    def plot_totaltorque(self):
        """
        Plot sum of the torque density
        """  
        ind = 24+self.nions-1
        tjxb = self.slices[0,0,:,7] #jxB torque from ini/end state
        tce = self.slices[0,0,:,ind] # collisional to el.
        tci1 = self.slices[0,0,:,ind+1] # collisional to ions

        if self.nions > 1:
            tci2 = self.slices[0,0,:,ind+2]
            if self.nions == 2:     
                plot_article(4,[self.rho, tjxb, tjxb+tce, tjxb+tci1, tjxb+tci2],['jxB','el.', 'i1', 'i2'],r'$\rho$', r'Torque density (N m^{-2})', self.infile_n)        
            else:
                tci3 = self.slices[0,0,:,ind+3]
                plot_article(5,[self.rho, tjxb, tjxb+tce, tjxb+tci1, tjxb+tci2, tjxb+tci3],\
                            ['jxB','el.', 'i1', 'i2', 'i3'],r'$\rho$', r'Torque density (N m^{-2})', self.infile_n)        
        else:
            plot_article(3,[self.rho, tjxb,  tjxb+tce, tjxb+tci1], ['jxB','el.', 'i1'], r'$\rho$', 'p (MW/$m^3$)', self.infile_n)


    def store_groupdis_to_ascii(self, fname):
        """
        Store total data distribution to ASCII
        """
        try:
            self.slices_summed.mean()
        except:
            self._group_beams()
        header=''

        for oo,eloo in enumerate(self.lab):
            header+=self.lab[oo]+self.uni[oo]+'  '
                
        np.savetxt(fname, self.slices_summed, fmt='%.8e', header=header)


    def _evaluate_shellVol(self, rho_new):
        """
        Evaluate the volumes at rho_new
        """
        rho_old = np.linspace(0, 1., len(self._volumes))
        param_V_rho = interpolate.interp1d(rho_old, self._volumes)
        self.volumes = param_V_rho(rho_new)
        
    def _evaluate_shellArea(self, rho_new):
        """
        Evaluate the areas at rho_new
        """
        rho_old = np.linspace(0, 1., len(self._areas))
        param_A_rho = interpolate.interp1d(rho_old, self._areas)
        self.areas = param_A_rho(rho_new)

    def _group_beams(self):
        """
        Group data distribution forgetting about beam index
	     """
        self.slices_summed = np.sum(self.slices, axis=0)


    def _calculate_scalar(self):
        """
        Calculate scalar quantities: total power to ions, to electrons,
        total current induced, total angular momentum
        """
        try:
            np.mean(self.slices_summed)
        except:
            self._group_beams()
            
        ind = self.name_dict['pel']-1
        if self.slices.shape[2] == 27:
            ind=ind-1
        
        self.I_tot = np.dot(self.slices_summed[0,:,12], self.areas)
        self.pe    = np.dot(self.slices_summed[0,:,ind], self.volumes)
        self.pi1   = np.dot(self.slices_summed[0,:,ind+1], self.volumes)
        self.tjxb  = np.dot(self.slices_summed[0,:,7], self.volumes) 
        self.tore  = np.dot(self.slices_summed[0,:,ind+2 + self.nions-1], self.volumes)
        self.tori1 = np.dot(self.slices_summed[0,:,ind+2 + self.nions], self.volumes)

        if self.nions > 1:
            self.pi2   = np.dot(self.slices_summed[0,:,ind+2], self.volumes)
#            self.tori2 =  np.dot(self.slices_summed[0,:,ind+3 + self.nions], self.volumes)
            if self.nions>2:
                self.pi3   = np.dot(self.slices_summed[0,:,ind+3], self.volumes)
#                self.tori3 = np.dot(self.slices_summed[0,:,ind+4 + self.nions], self.volumes)

    def print_scalars(self):
        """
        Print the scalars
        """
        try:
            self.I_tot
        except:
            self._calculate_scalar()

        pextra = 0
        textra = 0
        print(" ")
        print("Total current induced    ", self.I_tot*1e-3, " kA")
        print("Total power to electrons ", self.pe*1e-6, " MW")
        print("Total power to ions      ", self.pi1*1e-6, " MW")
        if self.nions > 1:
            print("Total power to ion  2 ", self.pi2*1e-6, " MW")
            pextra += self.pi2*1e-6
            if self.nions > 2:
                print("Total power to ion  3 ", self.pi3*1e-6, " MW")
                pextra += self.pi3*1e-6
        print("Total power delivered    ", (self.pi1+self.pe+pextra)*1e-6, " MW")
        
        print("JxB torque ", self.tjxb, " Nm")
        print("Torque to electrons", self.tore, " Nm")
        print("Torque to ions", self.tori1, " Nm")
        if self.nions > 1:
            print("Torque to ion 2 ", self.tori2*1e-6, " Nm")
            textra += self.tori2*1e-6
            if self.nions > 2:
                print("Torque to ion 3 ", self.tori3*1e-6, " Nm")
                textra += self.tori3*1e-6
        print("Torque delivered    ", (self.tori1+self.tore+textra)*1e-6, " Nm")
#            
            
    def calc_eff(self):                
        """
        Method to compute the current-drive efficiency from the beam:
        eta = R0*n_e*I_CD/P
        """                
        try:
            self.I_tot.mean()
        except:
            self._calculate_scalar()
        #Computing the power coupled to plasma: Pini-Pend+Pres
        pp = ascot_particles.h5_particles(self.infile_n); pp._power_coupled()
        P = pp.pcoup
        
        R0 = self.infile['misc/geomCentr_rz'][0]
        # average density in 10^20
        ne = self.infile['plasma/1d/ne'][:]
        rhone = self.infile['plasma/1d/rho'][:]
        param_ne = interpolate.interp1d(rhone, ne, kind='linear')
        vol = self.volumes
        ne_f = param_ne(np.linspace(0,1,len(vol)))        
        ne_avg = np.trapz(ne_f*vol)*1e-20
        ne_avg /= np.sum(self.volumes)
        #############################
        self.eff = R0*ne_avg*np.abs(self.I_tot)/P
        print("R0 ", R0, 'm')
        print("ne avg", ne_avg, '10e20 m^-3')
        print("Ip ", self.I_tot, 'A')
        print("P  ", P, 'W')
        print("CD efficiency: ", self.eff, " [10^20 A / (W m^2)]")


    def _slowingdown_3D(self):
        """
        3D calculations for collisionality (rho, Te, taus):
        (1) at each rho computes Ec
        (2) at each fraction of Te (2Te, 3Te, etc) computes the taus
        """         
        try:
            self.param_ec
        except:
            self._ecrit()
            
        # volume-averaged temperature
        te = self.infile['plasma/1d/te'][:]
        rhote = self.infile['plasma/1d/rho'][:]
        param_te = interpolate.interp1d(rhote, te, kind='linear')
        vol = self.volumes
        te_f = param_te(np.linspace(0,1,len(vol)))        
        te_avg = np.trapz(te_f*vol)
        te_avg /= np.sum(self.volumes)        

        te_frac = np.linspace(0,10, num=10)
        Ec = self.param_ec(self.rho)
        ts = self.param_ts(self.rho)
        taus = np.zeros((len(self.rho), len(te_frac)), dtype=float)
        E0 = 85e3*1.602e-19
        e_x = np.linspace(0, E0, num=100)
        for ind_te, i in enumerate(te_frac):
            for ind_ec,j in enumerate(Ec):
                j=j*1.602e-19
                coeff = np.trapz(e_x**0.5/(j**1.5+e_x**1.5),\
                                np.linspace(i*te_avg*1.602e-19, E0, num=100))      
                ts_t = ts[ind_ec]
                taus[ind_ec, ind_te] = ts_t*coeff
                
        self.taus=taus                
                
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
        self.ec_mean = np.trapz(Ec[rho<1]*self.volumes)/np.sum(self.volumes)
        #Spitzer slowing-down time
        ts = 6.28e14*A*te**1.5/(Z**2*ne*17.)
        self.param_ts = interpolate.interp1d(rho,ts)


class TCV_1d(distribution_1d):
    
    def __init__(self, infile_n):
        distribution_1d.__init__(self, infile_n)

    def TCV_calc_FIBP(self, plot_flag, *args):
        """
        Compute the FIBP in TCV case
        """
        volumes = self.volumes
        rho = np.linspace(0, 1, num=len(volumes), dtype=float)
        self.fibp        = np.zeros(len(rho),dtype=float)

        #rho_edg = rho+(rho[-1]-rho[-2])*0.5
        if len(args)==0:
            #origins  = self.infile['inistate/origin'].value
            #part_rho = self.infile['inistate/rho'].value
            weight_a = self.infile['inistate/weight'].value
        else:
            shot = args[0]
            run = args[1]
            new_fname = '/home/vallar/ASCOT/runs/TCV/'+"{:03d}".format(shot)+'/bbnbi_'+"{:03d}".format(shot)+"{:03d}".format(run)+'.h5'
            print("File opened: ", new_fname)
            bbfile = h5py.File(new_fname)
            #origins  = bbfile['inistate/origin'].value
            #part_rho = bbfile['inistate/rho'].value
            weight_a = bbfile['inistate/weight'].value
            
        for i,r in enumerate(rho):
            if r==rho[-1]:
                continue
     
            weight = np.sum(weight_a)

            self.fibp[i] = weight/volumes[i]

        if plot_flag == 1:
            plot_article(1,[rho, self.fibp], [''], r'$\rho$', r'Fast ion birth profile $1/(s\cdot m^3)$')

    
    def TCV_plot_all(self, *args):
        """
        Plot quantities without different ion species (i.e. j,p,ecc.)
        """
        if len(args)==0:
            y_ind = np.array([(self.name_dict[t]-1) for t in self.name_dict.keys()])
            n_row = 5
            n_col = 5
        else:
            y_ind = args[0]  
            if len(y_ind)<3:
                n_row = 1
            else:
                n_row = 3
            n_col = int(len(args[0])/n_row)
        x = self.rho
        y = np.zeros((np.size(y_ind), np.size(x)))
        y = self.slices[0,0,:,y_ind]
        n_el = np.size(y_ind)
        ylabels = np.array(self.lab[y_ind])
        yunits  = np.array(self.uni[y_ind])
        print(x.shape)
        print(y.shape)
        fig, ax = plt.subplots(n_row,n_col)
        for i in range(n_el):
            ind_row=i%n_row
            ind_col=i/n_col
            ax[ind_row,ind_col].plot(x, y[i,:])                    
            ax[ind_row,ind_col].set_xlabel('rho')
            ax[ind_row,ind_col].set_ylabel(yunits[i])
            ax[ind_row,ind_col].set_title(ylabels[i])

        plt.show()        


class SA_1d(distribution_1d):
    colours = ['g','c','b','k','r']

    def __init__(self, infile_n):
        distribution_1d.__init__(self, infile_n)
        
    def SA_plot_profs(self, *args):
        """
        Method to plot quantities without different ion species (i.e. j,p,ecc.)
        """
        if len(args)==0:
            y_ind = np.array([(self.name_dict[t]-1) for t in self.name_dict.keys()])
            n_row = 5
            n_col = 5
        else:
            y_ind = args[0]  
            if len(y_ind)<3:
                n_row = 1
            else:
                n_row = 3
            n_col = int(len(args[0])/n_row)
        x = self.rho
        y = np.zeros((np.size(y_ind), 26, np.size(x)))
        y = self.slices[:,0,:,y_ind]
        n_el = np.size(y_ind)
        ylabels = np.array(self.lab[y_ind])
        yunits  = np.array(self.uni[y_ind])

        
        fig, ax = plt.subplots(n_row,n_col)
        for i in range(n_el):
            ind_row=i%n_row
            ind_col=i/n_col
            for j in range(self.n_inj/2):
                if j!=12:
                    ax[ind_row,ind_col].plot(x, y[i,2*j,:]+y[i,2*j+1,:], color=colours[0])                    
                    #ax[ind_row,ind_col].plot(x, y[2*j,:,i]+y[2*j+1,:,i], color=colours[0])
                else:
                    ax[ind_row,ind_col].plot(x,y[i,2*j,:])
                    ax[ind_row,ind_col].plot(x,y[i,2*j+1,:])
                    #ax[ind_row,ind_col].plot(x,y[2*j,:,i])
                    #ax[ind_row,ind_col].plot(x,y[2*j+1,:,i])
            ax[ind_row,ind_col].set_xlabel('rho')
            ax[ind_row,ind_col].set_ylabel(yunits[i])
            ax[ind_row,ind_col].set_title(ylabels[i])

            
        plt.show()

    def _calc_originbeam(self):
        """
        Method to translate from origin in bbnbi.h5 file to beam identification
        Now only with JT60-SA (and still to recognize A and B units)
        """
        #self.beamorigindict={'1':[45,46]    , '2':[47,48]     ,'3':[133,134]   ,'4':[135,136],\
        #                     '5':[221,222]  , '6':[223,224]   ,'7':[309,310]   ,'8':[311,312],\
        #                     '9':[3637,3638], '10':[3639,3640],'13':[5253,5254],'14':[5255,5256],\
        #                     'NNBI_U':[3031], 'NNBI_L':[3032]}
        #self.origindict={'45':'1_1', '46':'1_2', '47':'2_1', '48':'2_2',\
        #                 '133':'3_1', '134':'3_2', '135':'4_1', '136':'4_2',\
        #                 '221':'5_1', '222':'5_2', '223':'6_1', '224':'6_2',\
        #                 '309':'7_1', '310':'7_2', '311':'8_1', '312':'8_2',\
        #                 '3637':'9_1', '3638':'9_2', '3639':'10_1', '3640':'10_2',\
        #                 '5253':'13_1', '5254':'13_2', '5255':'14_1', '5256':'14_2',\
        #                 '3031':'NNB_1', '3032':'NNB_2'}

#        self.beamlabel=['1','2','3','4','5','6','7','8','9','10',\
#                        '13','14','NNB_U','NNB_L']
#        self.beamlabel_full={'1':[45, 46],'2':[47, 48],\
#                            '3':[133, 134],'4':[135, 136],\
#                            '5':[221, 222],'6':[223, 224],\
#                            '7':[309, 310],'8':[311, 312],\
#                            '9':[3637, 3638],'10':[3639, 3640],\
#                            '13':[5253, 5254],'14':[5255, 5256],\
#                            '99':[3031],'101':[3032]}
#        self.beamorigin_full = np.array([45,46,47,48,133,134,135,136,\
#                           221,222,223,224,309,310,311,312,\
#                           3637,3638,3639,3640,5253,5254,5255,5256,\
#                           3031,3032], dtype=int)
        
        self.beamorigindict_full = {'45':1 ,  '46':1,    '47':2,'48':2,\
                                '133':3,  '134':3,\
                                '135':4,  '136':4,   '221':5,'222':5,\
                                '223':6,  '224':6,   '309':7,'310':7,\
                                '311':8,  '312':8,   '3637':9,'3638':9,\
                                '3639':10,'3640':10, \
                                '5253':13,'5254':13, '5255':14,'5256':14,\
                                '3031':99,'3032':101}
#        self.beamlabel_no78 = ['1','2','3','4','5','6','9','10',\
#                        '13','14','NNB_U','NNB_L']

        #orderedorigins is the array where you can retrieve the correct beam used in the distributions
        # so in orderedorigins we will find e.g. [13,7,9,11,15,4,8,9,...]
        # and you get data of beam 1 using orderedorigins[0], etc.
        # origindict is the dictionary where you can find in the keys the beam
        # number and the value is the element in orderedorigin
#        self.orderedorigins=np.zeros(len(self._h5origins), dtype=int)
#        self.origindict = {}
#        j=0
#        for kk in self._h5origins:
#            for jj, el in enumerate(self.beamorigin_full):
#                if kk==el:
#                    self.orderedorigins[j]=jj
#                    j=j+1         
#                    break
#            for i,el in enumerate(self.beamlabel_full.keys()):
#                if kk in self.beamlabel_full[el]:
#                    if el not in self.origindict.keys():
#                        self.origindict[el] = jj
#                    else:
#                        self.origindict[el] = np.append(self.origindict[el], jj)

    def SA_calculate_scalar(self):
        """
        Method to calculate scalar quantities: total power to ions, to electrons,
        total current induced, total angular momentum
        """
        
        self._evaluate_shellArea(self.rho)
        self._evaluate_shellVol (self.rho)
        self._calc_originbeam  ()
        self.I_tot = 0
        self.pe=0
        self.pi1 = 0
        self.tor = 0
        self.i_beams  = np.zeros(self.n_inj/2+1, dtype=float)
        self.pi_beams = np.zeros(self.n_inj/2+1, dtype=float)
        self.pe_beams = np.zeros(self.n_inj/2+1, dtype=float)
        self.tor_beams = np.zeros(self.n_inj/2+1, dtype=float)
        for jj in range(self.n_inj/2-1):
            ind1=self.orderedorigins[2*jj]
            ind2=self.orderedorigins[2*jj+1]
            # CURRENT
            tmp_j  = self.slices[ind1,0,:,2]+self.slices[ind2,0,:,2]
            #plt.plot(self.rho, tmp_j, self.rho, self.areas)
            #plt.show()
            tmp_I  = np.dot(tmp_j, self.areas)
            self.i_beams[jj] = tmp_I
            self.I_tot += tmp_I

            # POWER TO electrons
            tmp_Penorm = self.slices[ind1,0,:,21]+self.slices[ind2,0,:,21]
            tmp_Pe = np.dot(tmp_Penorm, self.volumes)
            self.pe_beams[jj] = tmp_Pe
            self.pe += tmp_Pe

            # POWER TO IONS 1
            tmp_Pinorm = self.slices[ind1,0,:,22]+self.slices[ind2,0,:,22]
            tmp_Pi = np.dot(tmp_Pinorm, self.volumes)
            self.pi_beams[jj] = tmp_Pi
            self.pi1 += tmp_Pi

            #TORQUE TOTAL
            tmp_tornorm = self.slices[ind1,0,:,5]+self.slices[ind1,0,:,24]+\
                          self.slices[ind2,0,:,5]+self.slices[ind2,0,:,24]
            tmp_tor = np.dot(tmp_tornorm, self.volumes)
            self.tor_beams[jj] = tmp_tor
            self.tor += tmp_tor
        #==============================================================
        # NNB
        ind1=self.orderedorigins[-2]
        ind2=self.orderedorigins[-1]
        # CURRENT
        tmp_j  = self.slices[ind1,0,:,2]
        #plt.plot(self.rho, tmp_j, self.rho, self.areas)
        #plt.show()
        tmp_I  = np.dot(tmp_j, self.areas)
        self.i_beams[-2] = tmp_I
        self.I_tot += tmp_I

        tmp_j  = self.slices[ind2,0,:,2]
        #plt.plot(self.rho, tmp_j, self.rho, self.areas)
        #plt.show()
        tmp_I  = np.dot(tmp_j, self.areas)
        self.i_beams[-1] = tmp_I
        self.I_tot += tmp_I
        
        # POWER TO electrons
        tmp_Penorm = self.slices[ind1,0,:,21]
        tmp_Pe = np.dot(tmp_Penorm, self.volumes)
        self.pe_beams[-2] = tmp_Pe
        self.pe += tmp_Pe

        tmp_Penorm = self.slices[ind2,0,:,21]
        tmp_Pe = np.dot(tmp_Penorm, self.volumes)
        self.pe_beams[-1] = tmp_Pe
        self.pe += tmp_Pe
        # POWER TO IONS 1
        tmp_Pinorm = self.slices[ind1,0,:,22]
        tmp_Pi = np.dot(tmp_Pinorm, self.volumes)
        self.pi_beams[-2] = tmp_Pi
        self.pi1 += tmp_Pi
        
        tmp_Pinorm = self.slices[ind2,0,:,22]
        tmp_Pi = np.dot(tmp_Pinorm, self.volumes)
        self.pi_beams[-1] = tmp_Pi
        self.pi1 += tmp_Pi

        # TORQUE
        tmp_tornorm = self.slices[ind1,0,:,5]+self.slices[ind1,0,:,24]
        tmp_tor = np.dot(tmp_tornorm, self.volumes)
        self.tor_beams[-2] = tmp_tor
        self.tor += tmp_tor

        tmp_tornorm = self.slices[ind2,0,:,5]+self.slices[ind2,0,:,24]
        tmp_tor = np.dot(tmp_tornorm, self.volumes)
        self.tor_beams[-1] = tmp_tor
        self.tor += tmp_tor       
        #==============================================================

    def print_scalars(self):
        """
        Method to print the scalars
        """
        try:
            np.mean(self.i_beams)
        except:
            self.SA_calculate_scalar()
                

        for ii, el in enumerate(self.i_beams):
            print("Current from beam ", self.beamlabel[ii], " is ", el*1e-3, " kA")
        print("")
        for ii, el in enumerate(self.pe_beams):
            print("Pe from beam ", self.beamlabel[ii], " is ", el*1e-6, " MW")
        print("")
        for ii, el in enumerate(self.pi_beams):
            print("Pi from beam ", self.beamlabel[ii], " is ", el*1e-6, " MW")
        print("")
        for ii, el in enumerate(self.tor_beams):
            print("Tor from beam ", self.beamlabel[ii], " is ", el, " Nm")
            
        print(" ")
        print("Total current induced    ", self.I_tot*1e-3, " kA")
        print("Total power to electrons ", self.pe*1e-6, " MW")
        print("Total power to ions      ", self.pi1*1e-6, " MW")
        print("Total power delivered    ", (self.pi1+self.pe)*1e-6, " MW")
        print("Total torque ", self.tor, " Nm")

    def group_beams(self):
        """
        Method to group data distribution for beam type (PPERP, PPAR, NNB)
        """
        self._calc_originbeam  ()
                
        self.data_PPERP = np.zeros((len(self.rho), len(self.lab)),dtype=float)
        self.data_PPAR  = np.zeros((len(self.rho), len(self.lab)),dtype=float)
        self.data_NNB   = np.zeros((len(self.rho), len(self.lab)),dtype=float)

        for ii, el in enumerate(self._h5origins):
        #for ii, el in enumerate(['NNB_L', 'NNB_U']):
            beam=self.beamorigindict_full[str(el)]
            if beam in [7,8,9,10]:
                self.data_PPAR += self.slices[ii,0,:,:]                  
            elif beam in [99, 101]:
                self.data_NNB += self.slices[ii, 0,:,:]
            else:
                self.data_PPERP += self.slices[ii,0,:,:]
 
    def plot_groupcurrent(self):
        """
        method to plot the data produced with group_beams 
        """
        try:
            self.data_PPAR.mean()
        except:
            self.group_beams()

        i_t_ppar = -1.*self.data_PPAR [:,12]*1e-3
        i_t_pper = -1.*self.data_PPERP[:,12]*1e-3
        i_t_nnb  = -1.*self.data_NNB  [:,12]*1e-3
        i_tot  = i_t_ppar+i_t_pper+i_t_nnb
        self.Itot_prof = i_tot
        
        n_ppar = self.data_PPAR [:,0]
        n_pper = self.data_PPERP[:,0]
        n_nnb  = self.data_NNB  [:,0]
        n_tot  = n_ppar + n_pper + n_nnb
        self.ntot_prof = n_tot

        pe_ppar = self.data_PPAR [:,21]*1e-3
        pe_pper = self.data_PPERP[:,21]*1e-3
        pe_nnb  = self.data_NNB  [:,21]*1e-3
        pe_tot  = pe_ppar+pe_pper+pe_nnb
        self.petot_prof = pe_tot
        
        pi_ppar = self.data_PPAR [:,22]*1e-3
        pi_pper = self.data_PPERP[:,22]*1e-3
        pi_nnb  = self.data_NNB  [:,22]*1e-3      
        pi_tot  = pi_ppar+pi_pper+pi_nnb
        self.pitot_prof = pi_tot

#        thn_ppar=self.data_PPAR  [:,18]
#        thn_pper=self.data_PPERP [:,18]
#        thn_nnb=self.data_NNB    [:,18]

        press_ppar = self.data_PPAR [:,13]*1e-3
        press_pper = self.data_PPERP[:,13]*1e-3
        press_nnb  = self.data_NNB  [:,13]*1e-3
        p_tot      = press_ppar+press_pper+press_nnb

        tor_i_ppar = self.data_PPAR [:,5]+self.data_PPAR [:,24]
        tor_i_pper = self.data_PPERP[:,5]+self.data_PPERP[:,24]
        tor_i_nnb  = self.data_NNB  [:,5]+self.data_NNB  [:,24]
#        tor_i_tot  = tor_i_ppar+tor_i_pper+tor_i_nnb
        
#        labels=['P-T','P-P','N-NB']
        labels2=['P-T','P-P','N-NB','TOT']
        plot_article(4,[self.rho, i_t_ppar, i_t_pper, i_t_nnb, i_tot],labels2,r'$\rho$', 'j (kA/$m^2$)', self.infile_n)
#        plot_article(3,[self.rho, n_ppar, n_pper, n_nnb],labels,r'$\rho$', 'n ($m^{-3}$)')
#        plot_article(4,[self.rho, pe_ppar, pe_pper, pe_nnb, pe_tot],labels2,r'$\rho$', '$P_e$ (kW/$m^3$)', self.infile_n)
#        plot_article(4,[self.rho, pi_ppar, pi_pper, pi_nnb, pi_tot],labels2,r'$\rho$', '$P_i$ (kW/$m^3$)', self.infile_n)
#        plot_article(4,[self.rho, press_ppar, press_pper, press_nnb, p_tot],labels2,r'$\rho$', '$p$ (kPa)', self.infile_n)
#        plot_article(4,[self.rho, tor_i_ppar, tor_i_pper, tor_i_nnb, tor_i_tot],labels2,r'$\rho$', 'Torque density (N $m^{-2}$)')

#        plot_article(3,[self.rho, tor_el_ppar, tor_el_pper, tor_el_nnb],labels,r'$\rho$', 'Torque to electrons (N $m^{-2}$)')
#        plot_article(4,[self.rho, i_t_ppar, i_t_pper, i_t_nnb, i_tot],labels2,r'$\rho$', 'Shielded current (kA/$m^{2}$)', self.infile_n)

        #plot_article(3,[self.rho, thn_ppar, thn_pper, thn_nnb],labels,r'$\rho$', 'Th. n (1/$m^3$)')

       
class distribution_2d:

    def __init__(self, infile_n):
        if os.path.isfile(infile_n) is False:
            print("File ", infile_n, " doesn't exists!")
            raise Exception()

        self.infile=h5py.File(infile_n, 'r')
        self.infile_n = infile_n
        if "rzPitchEdist" not in self.infile['distributions'].keys():
            print("No rzPitchE dist in ", self.infile_n)
        if "rhoPhiPEdist" not in self.infile['distributions'].keys():
            print("No rhoPhiPitchE dist in ", self.infile_n)
        if "rzMuEdist" not in self.infile['distributions'].keys():
            print("No rzMuE dist in ", self.infile_n)
        
        self.id = self.infile_n[-9:-3]
        self._readwall()
        
    def plot_space(self):
        try:
            self.f_Ep_int.mean()
        except:
            self._integrate_Ep()
       
        self.zplot = self.f_Ep_int
        if 'R' in self.dict_dim.keys() and 'z' in self.dict_dim.keys():
            self.xplot = self.dict_dim['R']
            self.yplot = self.dict_dim['z'] 
            self._plot_2d('R [m]', 'z [m]', wall=1, surf=1)
        elif 'rho' in self.dict_dim.keys() and 'phi' in self.dict_dim.keys():
            self.xplot = self.dict_dim['rho']
            self.yplot = self.dict_dim['phi'] 
            self._plot_2d('rho', 'phi', wall=0)
            

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
            
        int_E = np.trapz(dist_toint, self.dict_dim['E'], axis=0)
        self.f_Ep_int = np.trapz(int_E, self.dict_dim['pitch'], axis=0)


    def _integrate_spaceE(self):
        """
        Hidden method to integrate over (space,E)
        """
        self.f_spaceE_int = self._integrate_spacex('E', axis=0)


    def _integrate_spacep(self):
        """
        hidden method to integrate over (space,pitch)
        """
        self.f_spacep_int = self._integrate_spacex('pitch', axis=1)


    def _integrate_spacemu(self):
        """
        hidden method to integrate over (space,mu)
        """
        self.f_spacemu_int = self._integrate_spacex('mu')


    def _integrate_spacex(self, x, axis):
        """
        Hidden method to integrate over space and something else (pitch, E, mu...)
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        return np.trapz(self.f_space_int, self.dict_dim[x], axis=axis)        

    def plot_spacep(self):
        """
        plot 1D (energy, int_space (int_pitch (fdist)))
        """
        try:
            self.f_spacep_int.mean()
        except:
            self.integrate_spacep()
        
        self.xplot = self.dict_dim['E']/1.6e-19
        self.yplot = self.f_spacep_int
            
        self._plot_1d('E [keV]', "Normalized f")


    def plot_spaceE(self):
        """
        plot 1D (pitch, int_space (int_E (fdist)))
        """
        try:
            self.f_spaceE_int.mean()
        except:
            self._integrate_spaceE()
        
        self.xplot = self.dict_dim['pitch']
        self.yplot = self.f_spaceE_int

        self._plot_1d(r'$\xi$ ($\frac{v_\parallel}{v}$)', "Normalized f")


    def plot_Epitch(self):
        """
        plot 2D (pitch, energy, int_space(fdist))
        """
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()

        self.xplot = self.dict_dim['pitch']
        self.yplot = self.dict_dim['E']/1.6e-19
        self.zplot = self.f_space_int
        self._plot_2d('pitch', 'E')

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
        self._plot_2d('mu', 'E', wall=0)
        
    def write_pitchE(self):
        try:
            self.f_space_int.mean()
        except:
            self._integrate_space()
        self._write('pitch','E', self.f_space_int, units=['adimensional', 'J'])

    def _plot_1d(self, xlab, ylab):
        """
        Hidden method to plot 1D functions
        """

        fig = plt.figure()
        fig.suptitle(self.id)
        ax  = fig.add_subplot(111)
        ax.plot(self.xplot, self.yplot, 'k', linewidth=2.3)

        ax.set_xlabel(xlab), ax.set_ylabel(ylab)
        
        plt.show()

    def _plot_2d(self, xlab, ylab, **kwargs):
        """
        Hidden method to plot the 2D distribution
        wall: set to 1 if wall needed to plot (i.e. RZ function)
        """

        flag_dict = kwargs
        title = 'Normalized f'
        fig = plt.figure()
        tit=self.id
        ax  = fig.add_subplot(111)
        CS  = ax.contourf(self.xplot, self.yplot, self.zplot, 50, cmap=my_cmap)
        plt.colorbar(CS)
        #ax.plot([np.min(self.xplot),np.max(self.xplot)], [5e5, 5e5], linewidth=3., color='k')
               
        if 'wall' in flag_dict.keys() and flag_dict['wall']==1:
            ax.plot(self.R_w, self.z_w, 'k', linewidth=2)
        if 'surf' in flag_dict.keys() and flag_dict['surf']==1:
            self._plot_RZsurf(ax)            
        if 'title' in flag_dict.keys():
            tit+= ' '+flag_dict['title']
        fig.suptitle(tit)
        
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.grid('on')
        fig.tight_layout()
        plt.show()
        if 'fname' in flag_dict.keys():
            plt.savefig(flag_dict['fname'], bbox_inches='tight')
            
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

        self.info = '' + self.infile_n + ' ' + ' matteo.vallar@igi.cnr.it ' + \
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


    def _computenorm(self):
        """
        calculates the norm of the function
        """
        if "pitch" in self.dict_dim.keys():
            try:
                self.f_spacep_int()
            except:
                self._integrate_spacep()
                
            self.norm = np.trapz(self.f_spacep_int, self.dict_dim['E'])
            #print "NORM = ", self.norm
        
        elif "mu" in self.dict_dim.keys():
            try:
                self.f_spacemu_int()
            except:
                self._integrate_spacemu()
                
            self.norm = np.trapz(self.f_spacemu_int, self.dict_dim['E'])
            #print "NORM = ", self.norm        

    def build_fdist(self):
        """
        Method to read the ordinates of the 2d distribution
        """
        try:
            self.dict_dim.keys()
        except:
            print("No dictionary of dimensions created")
            raise ValueError

        self.norm = 1
        # 6th dimension is the one labelling the beams
        tmp = self.dist_h5['ordinate'].value
        fdist = np.sum(tmp, axis=0)[:,:,:,:,:,0] #time, E,pitch,z,r
        self.fdist_notnorm = fdist
        self._computenorm()
        #self.fdist_norm = self.fdist_notnorm/self.norm

    def collect_dim(self):
        """
        methods to read the dictionary of the dimensions of the 4D distributions and store it
        """
        self._read_dim()
        self._fix_dim()

    def _read_dim(self):
        """
        Hidden method to read the abscissae
        """
        self.shape_dim = self.dict_dim.copy()
        for i, dim in enumerate(self.dict_dim.keys()):
            self.dict_dim[dim] = self.dist_h5['abscissae/dim'+str(i+1)].value
            self.shape_dim[dim] = np.size(self.dict_dim[dim])-1

    def _fix_dim(self):
        """
        Hidden method to make the abscissae the correct length (same as ordinate)
        """
        try:
            self.dict_dim.keys()
        except:
            self._read_dim()

        for dim in self.dict_dim.keys():
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

    def _plot_RZsurf(self, ax):
        try:
            self.RZsurf.mean()
        except:
            self._RZsurf()
            
        CS = ax.contour(self.Rsurf, self.zsurf, self.RZsurf, [0.2, 0.4, 0.6, 0.8, 1.0], colors='k')
        plt.clabel(CS, inline=1, fontsize=10)   
        
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
        
    def _integrate_space_enslice(self, sliceind):
        """
        Function to integrate over (R,z) on a E defined
        """
        dist_toint = self.fdist_notnorm[0,sliceind,:,:,:]/self.norm

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
        makes a plot of E, pitch on a R,z position
        """
        print(kwargs)
        ind_R = np.argmin(self.dict_dim['R']-sliceR < 0)
        ind_z = np.argmin(self.dict_dim['z']-slicez < 0)
        dist_toplot = self.fdist_notnorm[0,:,:,ind_R, ind_z]/self.norm

        self.xplot = self.dict_dim['pitch']
        self.yplot = self.dict_dim['E']/1.6e-19
        self.zplot = dist_toplot
        if 'fname' in kwargs.keys():
            self._plot_2d('pitch', 'E', \
                          title='R='+str(sliceR)+' z='+str(slicez), \
                          fname=kwargs['fname'])
        #else:
        #    self._plot_2d('pitch', 'E', title='R='+str(sliceR)+' z='+str(slicez))
        

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
        int_rhophi = int_rho[:,:,0]
        self.f_space_int = int_rhophi #E,pitch  

class frzmue(distribution_2d):
    
    def __init__(self, infile_n):
        """
        Module to initialise the distributions RZMUEDIST
        self.fdist['ordinate'] has the following shape: (ind beam, time, energy, pitch, z, R, #ions)
        """
        distribution_2d.__init__(self, infile_n)
        try:
            self.dist_h5 = self.infile['distributions/rzMuEdist'] 
        except:
            raise ValueError
        self.__name__ = 'frzpe'

        self.dict_dim = collections.OrderedDict([('R',[]),('z',[]),('mu',[]),('E',[]),('t',[])])
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


def plot_article(n_lines, data, data_labels, xlabel, ylabel, title):
        #=====================================================================================
        # SET TEXT FONT AND SIZE
        #=====================================================================================
        plt.rc('font', family='serif', serif='Palatino')
        #plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=====================================================================================
        col=['k','r','g','b']
        fig=plt.figure()
        ax=fig.add_subplot(111)
        for i in range(n_lines):
            ax.plot(data[0], data[i+1], label=str(data_labels[i]), linewidth=3, color=col[i])
        ax.plot(data[0], np.zeros(len(data[0])), 'k',linewidth=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        #ax.set_ylim([0,250])
        #ax.plot([0.85,0.85],[min(ax.get_ybound()), max(ax.get_ybound())],'k--', linewidth=3.)
        #=====================================================================================
        # ADJUST SUBPLOT IN FRAME
        #=====================================================================================
        plt.subplots_adjust(top=0.95,bottom=0.12,left=0.15,right=0.95)
        #=====================================================================================
        plt.show()
        #=====================================================================================
        # SET TICK LOCATION
        #=====================================================================================

        # Create your ticker object with M ticks
        M = 4
        yticks = ticker.MaxNLocator(M)
        xticks = ticker.MaxNLocator(M)
        # Set the yaxis major locator using your ticker object. You can also choose the minor
        # tick positions with set_minor_locator.
        ax.yaxis.set_major_locator(yticks)
        #ax.yaxis.set_minor_locator(yticks_m)
        ax.xaxis.set_major_locator(xticks)
        #=====================================================================================
        #ax.set_ylim([0,150])
        if data_labels[0]!='':
            ax.legend(loc='best')
        fig.tight_layout()

