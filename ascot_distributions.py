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
from matplotlib import ticker
from scipy import interpolate
import os.path, math
import collections

colours = ['k','g','c','b','r']

class distribution_1d:
    """
    Class for handling the <distributions> data

    METHODS:
    __init__(self, infile): getting basic variables, checking if rhoDist in hdf5 file
    rhodists(self): Method to get the data from the ascot file
    TCV_calc_FIBP(self, plot_flag, *args): Compute the FIBP in TCV case. Plot_flag=1 to plot, if given as *args shot and run it plots the FIBP from a different BBNBI h5 file
    TCV_plot_all(self): Plot quantities without different ion species (i.e. j,p,ecc.)


    plot_current(self): Plot the induced beam current density
    plot_power(self): Plot just the deposited power density
    plot_torque(self): Plot just the deposited torque density
    plot_totalcurrent(self): Plot sum of the induced beam current density
    plot_totalpower(self): Plot sum of the deposited power density
    plot_totaltorque(self): Plot sum of the torque density

    store_groupdis_to_ascii(self, fname): Store total data distribution to ASCII

    print_scalars(self): Print the scalars

    
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
        try:
            rho_new = self.infile['distributions/rhoDist/abscissae/dim1'].value
            self.volumes = self.infile['distributions/rhoDist/shellVolume'].value
            self.areas   = self.infile['distributions/rhoDist/shellArea'].value
            #self._evaluate_shellVol(rho_new)
            #self._evaluate_shellArea(rho_new)
        except:
            print "No /distributions/rhoDist/ in ", infile_n    
        
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
                          'th_n':18, 'th_e_n':19, 'th_torque':20, 'abs_ICRH':21, 'J.B':22 ,\
                          'pel':23, 'pi1':24,\
                          'ctor_el':25, 'ctor_i1':26} # The last two lines are affected by the number of ion species
        
        self.abscissae = {}
        self.abscissae = self.abscissae.fromkeys(self.infile['/distributions/rhoDist/abscissae'].keys(),0)
        for key in self.abscissae.keys():
            self.abscissae[key]=self.infile[tree_path+'/abscissae/'+str(key)].value
        self.rho = np.linspace(0,1,len(self.abscissae['dim1'])-1)
            
        ordinate = self.infile['/distributions/rhoDist/ordinate'].value
        
        #ADDING DIFFERENT ION SPECIES
        self.nions = 1
        if ordinate.shape[-1]>25:
            diff=ordinate.shape[-1]-25
            n_ions_more=diff/2
            self.nions += n_ions_more
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
        plot_article(1,[self.rho, i_tot],[''],r'$\rho$', 'j (kA/$m^2$)')

    def plot_totalpower(self):
        """
        Plot sum of the deposited power density
        """
        ind = 23-1
        pe = self.slices[0,0,:,ind]*1e-6
        pi1 = self.slices[0,0,:,ind+1]*1e-6
        if self.nions > 1:
            pi2 = self.slices[0,0,:,ind+2]*1e-6
            if self.nions == 2:
                plot_article(3,[self.rho, pe, pi1, pi2],['el.', 'i1', 'i2'],r'$\rho$', 'p (MW/$m^3$)')        
                
            elif self.nions == 3:
                pi3 = self.slices[0,0,:,ind+3]*1e-6
                plot_article(4,[self.rho, pe, pi1, pi2, pi3],['el.', 'i1', 'i2', 'i3'],r'$\rho$', 'p (MW/$m^3$)')        

        else:
            plot_article(2,[self.rho, pe, pi1],['el.', 'i1'], r'$\rho$', 'p (MW/$m^3$)')        

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
                plot_article(4,[self.rho, tjxb, tjxb+tce, tjxb+tci1, tjxb+tci2],['jxB','el.', 'i1', 'i2'],r'$\rho$', r'Torque density (N m^{-2})')        
            else:
                tci3 = self.slices[0,0,:,ind+3]
                plot_article(5,[self.rho, tjxb, tjxb+tce, tjxb+tci1, tjxb+tci2, tjxb+tci3],['jxB','el.', 'i1', 'i2', 'i3'],r'$\rho$', r'Torque density (N m^{-2})')        
        else:
            plot_article(3,[self.rho, tjxb,  tjxb+tce, tjxb+tci1], ['jxB','el.', 'i1'], r'$\rho$', 'p (MW/$m^3$)')


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


    def TCV_calc_FIBP(self, plot_flag, *args):
        """
        Compute the FIBP in TCV case
        """
        volumes = self.volumes
        rho = np.linspace(0, 1, num=len(volumes), dtype=float)
        self.fibp        = np.zeros(len(rho),dtype=float)

        rho_edg = rho+(rho[-1]-rho[-2])*0.5
        if len(args)==0:
            origins  = self.infile['inistate/origin'].value
            part_rho = self.infile['inistate/rho'].value
            weight_a = self.infile['inistate/weight'].value
        else:
            shot = args[0]
            run = args[1]
            new_fname = '/home/vallar/ASCOT/runs/TCV/'+"{:03d}".format(shot)+'/bbnbi_'+"{:03d}".format(shot)+"{:03d}".format(run)+'.h5'
            print "File opened: ", new_fname
            bbfile = h5py.File(new_fname)
            origins  = bbfile['inistate/origin'].value
            part_rho = bbfile['inistate/rho'].value
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

        self.I_tot = np.dot(self.slices_summed[0,:,12], self.areas)
        self.pe    = np.dot(self.slices_summed[0,:,22], self.volumes)
        self.pi1   = np.dot(self.slices_summed[0,:,23], self.volumes)
        self.tjxb  = np.dot(self.slices_summed[0,:,7], self.volumes) 
        self.tore  = np.dot(self.slices_summed[0,:,24 + self.nions-1], self.volumes)
        self.tori1 = np.dot(self.slices_summed[0,:,24 + self.nions], self.volumes)

        if self.nions > 1:
            self.pi2   = np.dot(self.slices_summed[0,:,24], self.volumes)
            self.tori2 =  np.dot(self.slices_summed[0,:,24 + self.nions+1], self.volumes)
            if self.nions>2:
                self.pi3   = np.dot(self.slices_summed[0,:,25], self.volumes)
                self.tori3 = np.dot(self.slices_summed[0,:,24 + self.nions+2], self.volumes)


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
        print " "
        print "Total current induced    ", self.I_tot*1e-3, " kA"
        print "Total power to electrons ", self.pe*1e-6, " MW"
        print "Total power to ions      ", self.pi1*1e-6, " MW"
        if self.nions > 1:
            print "Total power to ion  2 ", self.pi2*1e-6, " MW"
            pextra += self.pi2*1e-6
            if self.nions > 2:
                print "Total power to ion  3 ", self.pi3*1e-6, " MW"
                pextra += self.pi3*1e-6
        print "Total power delivered    ", (self.pi1+self.pe+pextra)*1e-6, " MW"
        
        print "JxB torque ", self.tjxb, " Nm"
        print "Torque to electrons", self.tore, " Nm"
        print "Torque to ions", self.tori1, " Nm"
        if self.nions > 1:
            print "Torque to ion 2 ", self.tori2*1e-6, " Nm"
            textra += self.tori2*1e-6
            if self.nions > 2:
                print "Torque to ion 3 ", self.tori3*1e-6, " Nm"
                textra += self.tori3*1e-6
        print "Torque delivered    ", (self.tori1+self.tore+textra)*1e-6, " Nm"
                

class distribution_2d:

    def __init__(self, infile_n):
        if os.path.isfile(infile_n) is False:
            print "File ", infile_n, " doesn't exists!"
            raise Exception()

        self.infile=h5py.File(infile_n, 'r')
        self.infile_n = infile_n
        if "rzPitchEdist" not in self.infile['distributions'].keys():
            print "No rzPitchE dist in ", self.infile_n
        if "rhoPhiPEdist" not in self.infile['distributions'].keys():
            print "No rhoPhiPitchE dist in ", self.infile_n

    def plot_RZ(self):
        self.xplot = self.dict_dim['R']
        self.yplot = self.dict_dim['z']
        int_E = np.trapz(self.fdist[0,:,:,:,:], self.dict_dim['E'], axis = 0)
        int_pitchE = np.trapz(int_E, self.dict_dim['pitch'], axis=0)
        self.zplot = int_pitchE
        self._plot_2d('R', 'z', wall=1, norm=1)

    def integrate_RZ_rhophi(self):
        """
        Function to integrate over (R,z) or (rho,phi)
        """
        if 'R' in self.dict_dim.keys():
            multiply_R = 2*math.pi*np.tile(self.dict_dim['R'], (20,10,20,1))
            dist_toint = multiply_R*self.fdist[0,:,:,:,:]
            int_R   = np.trapz(dist_toint, self.dict_dim['R'], axis = -1)
            int_Rz  = np.trapz(int_R     , self.dict_dim['z'], axis = -1)
            self.f_RZ_int = int_Rz #E,pitch
        elif 'rho' in self.dict_dim.keys():
            int_rho  = np.sum(self.fdist[0,:,:,:,:], axis = -1)
            int_rhophi  = np.sum(int_rho, axis = -1)
            self.f_rhophi_int = np.transpose(int_rhophi) #E,pitch        

    def integrate_RZE(self):
        """
        Function to integrate over (R,z,E) or (rho,phi,E)
        """
        try:
            if 'R' in self.dict_dim.keys():
                self.f_RZ_int.mean()
            else:
                self.f_rhophi_int.mean()
        except:
            self.integrate_RZ_rhophi()
            
        if 'R' in self.dict_dim.keys():
            self.f_RZE_int = np.trapz(self.f_RZ_int, self.dict_dim['E'], axis=0)
            self.f_RZE_int = self.f_RZE_int
        else:
            self.rhophiE_int = np.sum(self.f_rhophi_int, axis=1)

    def integrate_RZp(self):
        """
        Function to integrate over (R,z,p) or (rho,phi,p)
        """
        try:
            if 'R' in self.dict_dim.keys():
                self.f_RZ_int.mean()
            else:
                self.f_rhophi_int.mean()
        except:
            self.integrate_RZ_rhophi()
            
        if 'R' in self.dict_dim.keys():
            self.f_RZp_int = np.trapz(self.f_RZ_int, self.dict_dim['pitch'], axis=1)
            self.f_RZp_int = self.f_RZp_int
        else:
            self.rhophip_int = np.sum(self.f_rhophi_int, axis=0)
            el_vol = self.dvol['rho']*self.dvol['phi']*self.dvol['pitch']
            self.f_rhophip_int = self.f_rhophip_int*el_vol/self.vtot

        

    def plot_RZp(self, norm):
        try:
            if 'R' in self.dict_dim.keys():
                self.f_RZp_int.mean()
            else:
                self.rhophip_int.mean()
        except:
            self.integrate_RZp()
        
        self.xplot = self.dict_dim['E']/1.6e-19
        if 'R' in self.dict_dim.keys():
            self.yplot = self.f_RZp_int
        else:
            self.yplot = self.f_rhophip_int
            
        self._plot_1d('Energy', norm=norm)



    def plot_RZE(self, norm):
        try:
            if 'R' in self.dict_dim.keys():
                self.f_RZE_int.mean()
            else:
                self.rhophiE_int.mean()
        except:
            self.integrate_RZE()
        
        self.xplot = self.dict_dim['pitch']
        if 'R' in self.dict_dim.keys():
            self.yplot = self.f_RZE_int
        else:
            self.yplot = self.f_rhophiE_int
            
        self._plot_1d('pitch', norm=norm)


    def plot_Epitch(self):
        try:
            if 'R' in self.dict_dim.keys():
                self.f_RZ_int.mean()
            else:
                self.f_rhophi_int.mean()
        except:
            self.integrate_RZ_rhophi()


        self.xplot = self.dict_dim['pitch']
        self.yplot = self.dict_dim['E']/1.6e-19

        if 'R' in self.dict_dim.keys():
            self.zplot = self.f_RZ_int
        else:
            self.zplot = self.f_rhophi_int
        self._plot_2d('pitch', 'E', wall=0, norm=1)


    def write_pitchE(self):
        try:
            self.f_RZ_int.mean()
        except:
            self.integrate_RZ_rhophi()
        self._write('pitch','E', self.f_RZ_int)

    def _plot_1d(self, xlab, norm):
        
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        title = 'Distribution '
        if norm==1:
            self.yplot = self.yplot/self.norm
            title += ' NORM'
        ax.plot(self.xplot, self.yplot, linewidth=2.3)

        ax.set_xlabel(xlab)
        ax.set_title(title)
        
        plt.show()

    def _plot_2d(self, xlab, ylab, wall, norm):
        flag_dict = {'wall': wall, 'norm': norm}
        title = 'Distribution function'
        
        if flag_dict['norm']==1:
            self.zplot = self.zplot/self.norm
            title+= ' NORMALIZED'
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        CS  = ax.contourf(self.xplot, self.yplot, self.zplot,50)
        CB  = plt.colorbar(CS)
        
        if os.path.isfile('input.wall_2d') & flag_dict['wall']!=0:
            data_w=np.loadtxt('input.wall_2d', skiprows=1)
            ax.plot(data_w[:,0], data_w[:,1], 'k', linewidth=2)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        plt.show()
    

    def _write(self, *args):
        """
        Method to write the distribution to dat file
        """
        x_labels = [args[i] for i in range(len(args)-1)]
        self.y = args[-1]
        fname = args[0]+args[1]+'.dat'
        self.header = ''
        for i in range(len(args)-1):
            self.header += args[i]+' '
        for i in range(len(args)-1):
            self.header += ' '+str(len(self.dict_dim[x_labels[i]]))

        with open(fname,'a') as f_handle:
            f_handle.write('# '+self.header+'\n')
            for lab in x_labels:
                np.savetxt(f_handle, self.dict_dim[lab].flatten())
            np.savetxt(f_handle, self.y.flatten())  

    def build_fdist(self):
        """
        Method to read the ordinates of the 2d distribution
        """
        try:
            self.dict_dim.keys()
        except:
            print "No dictionary of dimensions created"
            raise ValueError
            
        #self.fdist = np.zeros(1)
        #self.fdist = np.resize(self.fdist, [i for i in [len(k) for k in self.dict_dim.values()][::-1]])

        # 6th dimension is the one labelling the beams
        tmp = self.dist_h5['ordinate'].value
        self.fdist = np.sum(tmp, axis=0)[:,:,:,:,:,0] #time, E,pitch,z,r
        self.fdist_notnorm = self.fdist
        #computing normalisation of the function, integrating over all the dimensions
        self.integrate_RZp()
        self.norm = np.trapz(self.f_RZp_int, self.dict_dim['E'])
        self.fdist_norm = self.fdist_notnorm/self.norm


    def collect_dim(self):
        self._read_dim()
        self._fix_dim()

    def _read_dim(self):
        """
        Hidden method to read the abscissae
        """
        for i, dim in enumerate(self.dict_dim.iterkeys()):
            self.dict_dim[dim] = self.dist_h5['abscissae/dim'+str(i+1)].value

    def _fix_dim(self):
        """
        Hidden method to make the abscissae the correct length (same as ordinate)
        """
        try:
            self.dict_dim.keys()
        except:
            self._read_dim()
        tmp_dict = self.dict_dim
        for dim in self.dict_dim.keys():
            tmp_dim = tmp_dict[dim]
            self.dict_dim[dim] = np.linspace(min(tmp_dim), max(tmp_dim), len(tmp_dim)-1)


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
        
        self.dict_dim = collections.OrderedDict([('R',[]),('z',[]),('pitch',[]),('E',[]),('t',[])])

        self.collect_dim()
        self.build_fdist()


class frhophipe(distribution_2d):
    
    def __init__(self, infile_n):
        """
        Module to initialise the distributions rhoPhiPEdist
        self.fdist['ordinate'] has the following shape: (ind beam, time, energy, pitch, phi, rho, #ions)
        """
        distribution_2d.__init__(self, infile_n)
        try:
            self.dist_h5 = self.infile['distributions/rhoPhiPEdist']
        except:
            raise ValueError
        
        self.dict_dim = collections.OrderedDict([('rho',[]),('phi',[]),('pitch',[]),('E',[]),('t',[])])
        
        self.collect_dim()
        self.build_fdist()


def plot_article(n_lines, data, data_labels, xlabel, ylabel):
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
        #ax.set_ylim([0,240])
        #ax.plot([0.85,0.85],[min(ax.get_ybound()), max(ax.get_ybound())],'k--', linewidth=3.)
        #=====================================================================================
        # ADJUST SUBPLOT IN FRAME
        #=====================================================================================
        plt.subplots_adjust(top=0.95,bottom=0.12,left=0.15,right=0.95)
        #=====================================================================================

        #=====================================================================================
        # SET TICK LOCATION
        #=====================================================================================

        # Create your ticker object with M ticks
        M = 4
        yticks = ticker.MaxNLocator(M)
        yticks_m=ticker.MaxNLocator(M*2)
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
