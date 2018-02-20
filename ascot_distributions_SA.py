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
    __init__(self, infile) to store variables
    plot_profs(self,inlist) to plot values without different effect for different ions
        e.g. J, P, ...
    SA_plot_profs(): specific plots for JT60SA
    checkplot(self)        to plot values and check them
    calculate_scalar(self): Method to calculate scalar quantities: total power to ions, to electrons,
        total current induced, total angular momentum
    print_scalars(self): method to print the scalars
    store_scal_to_ascii(self): 	Method to store data scalars (all grouped in a single file) to ascii

    store_dis_to_ascii(self):	Method to store data distribution (one for each beam) to ASCII

    group_beams(self): 	Method to group data distribution for beam type (PPERP, PPAR, NNB)

    store_groupdis_to_ascii(self): 	Method to store data distribution (one for each BEAM TYPE) to ASCII

    plot_groupcurrent(self):         method to plot the data produced with group_beams 



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
        self.infile=h5py.File(self.infile)


    """
    def __init__(self, infile_n):
        self.infile=h5py.File(infile_n)
        self.infile_n = infile_n
        self.rho = self.infile['distributions/rhoDist/abscissae/dim1'].value
        try:
            self._volumes = self.infile['distributions/rhoDist/shellVolume'].value
            self._areas   = self.infile['distributions/rhoDist/shellArea'].value
            self._evaluate_shellVol()
            self._evaluate_shellArea()
        except:
            print "No /distributions/rhoDist/ in ", infile_n    

        self.rhodists()
        self.fibp_particles = 0
    
    def _print_options(self):
        """
        Hidden method to print the options of the ascot file
        """
        print self.infile['inputfiles/input.options'].value

    def _evaluate_shellVol(self):
        """
        Function to evaluate the volumes at the correct rho positions
        """
        rho_new = self.rho
        rho_old = np.linspace(0, 1., len(self._volumes))
        param_V_rho = interpolate.interp1d(rho_old, self._volumes)
        self.volumes = param_V_rho(rho_new)
        
    def _evaluate_shellArea(self):
        """
        Function to evaluate the Area at the correct rho positions
        """
        rho_new = self.rho
        rho_old = np.linspace(0, 1., len(self._areas))
        param_A_rho = interpolate.interp1d(rho_old, self._areas)
        self.areas = param_A_rho(rho_new)

        
    def plot_1d(self, n_lines, data, data_labels, xlabel, ylabel):
        #=====================================================================================
        # SET TEXT FONT AND SIZE
        #=====================================================================================
        plt.rc('font', family='serif', serif='Palatino')
        #plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=====================================================================================
        col=['r','g','b','k']
        fig=plt.figure()
        fig.suptitle(self.infile_n)
        ax=fig.add_subplot(111)
        for i in range(n_lines):
            ax.plot(data[0], data[i+1], label=str(data_labels[i]), linewidth=3, color=col[i])
        #ax.plot(data[0], np.zeros(len(data[0])), 'k',linewidth=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_ylim([0,2.5e19])
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
        ax.legend(loc='best')


    def SA_calc_FIBP(self, plot_flag, *args):
        volumes = self.volumes
        rho = np.linspace(0, 1, num=len(volumes), dtype=float)
        self.fibp        = np.zeros(len(rho),dtype=float)
        self.fibp_PNB    = np.zeros(len(rho),dtype=float)
        self.fibp_P_tang = np.zeros(len(rho),dtype=float)
        self.fibp_P_perp = np.zeros(len(rho),dtype=float)
        self.fibp_NNB    = np.zeros(len(rho),dtype=float)
        rho_edg = rho+(rho[-1]-rho[-2])*0.5
        if len(args)==0:
            origins  = self.infile['inistate/origin'].value
            part_rho = self.infile['inistate/rho'].value
            weight_a = self.infile['inistate/weight'].value
        else:
            shot = args[0]
            run = args[1]
            print shot, run
            new_fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/bbnbi_'+"{:03d}".format(shot)+"{:03d}".format(run)+'.h5'
            print "File opened: ", new_fname
            bbfile = h5py.File(new_fname)
            origins  = bbfile['inistate/origin'].value
            part_rho = bbfile['inistate/rho'].value
            weight_a = bbfile['inistate/weight'].value
            
        for i,r in enumerate(rho):
            if r==rho[-1]:
                continue
            ind = [(part_rho>rho_edg[i]) &  (part_rho<rho_edg[i+1])]
            ind_P = [(part_rho>rho_edg[i]) & (part_rho<rho_edg[i+1]) & (origins != 3031) & (origins != 3032)]
            ind_P_perp = [(part_rho>rho_edg[i]) & (part_rho<rho_edg[i+1]) & \
                          (origins != 3031) & (origins != 3032) & \
                          (origins != 309) & (origins != 310) & \
                          (origins != 312) & (origins != 311) & \
                          (origins != 3637) & (origins != 3638) & \
                          (origins != 3639) & (origins != 3640)\
                          ]
        
            weight = np.sum(weight_a[ind])
            weight_P = np.sum(weight_a[ind_P])
            weight_PP = np.sum(weight_a[ind_P_perp])

            self.fibp[i] = weight/volumes[i]
            self.fibp_PNB[i]=weight_P/volumes[i]
            self.fibp_P_perp[i]=weight_PP/volumes[i]
            self.fibp_P_tang[i]=self.fibp_PNB[i]-self.fibp_P_perp[i]
            self.fibp_NNB[i]=self.fibp[i]-self.fibp_PNB[i]

        self.fibp_particles = np.dot(self.fibp, volumes)
        if plot_flag == 1:
            self.plot_1d(4, [rho, self.fibp, self.fibp_P_perp, self.fibp_P_tang, self.fibp_NNB], \
                     ['TOT','P-perp','P-tang','N'], r'$\rho$', r'Fast ion birth profile $1/(s\cdot m^3)$')

    def TCV_calc_FIBP(self, plot_flag, *args):
        volumes = self.volumes
        rho = np.linspace(0, 1, num=len(volumes), dtype=float)
        self.fibp        = np.zeros(len(rho),dtype=float)
        self.fibp_PNB    = np.zeros(len(rho),dtype=float)
        self.fibp_P_tang = np.zeros(len(rho),dtype=float)
        self.fibp_P_perp = np.zeros(len(rho),dtype=float)
        self.fibp_NNB    = np.zeros(len(rho),dtype=float)
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
            self.plot_1d(1, [rho, self.fibp], ['TOT'], \
                         r'$\rho$', r'Fast ion birth profile $1/(s\cdot m^3)$')
    
    def rhodists(self):
            
        tree_path = '/distributions/rhoDist/'

        self.name_dict = {'n':1, 'e_den':2,\
                          'jpar':3, 'jperp':4, \
                          'jxB':5, 'jxBstate':6, \
                          'CX_ionsource':7, 'CX_ionensource':8,\
                          'CX_neutsource':9, 'CX_neutensource':10,\
                          'torque':11,'par_e_den':12,\
                          'tot_tor_j':13, 'ptot':14, 'ppar':15, 'pperp':16,\
                          'flr_torque':17,\
                          'th_n':18, 'th_e_n':19, 'th_torque':20, 'abs_ICRH':21,\
                          'j.b':22,\
                          'pel':23, 'pi1':24,\
                          'ctor_el':25, 'ctor_i1':26\
                          }

        self.thlist   = ['n', 'e_den', 'th_n', 'th_e_n', 'th_torque' ]
        self.plist    = ['ptot', 'ppar', 'pperp']
        self.jlist    = ['jpar', 'jperp','jxB', 'tot_tor_j']
        self.pdeplist = ['pel', 'pi1', 'pi2', 'pi3']
        self.torlist  = ['ctor_el', 'ctor_i1', 'ctor_i2', 'ctor_i3'] 


        self.abscissae = {}
        # these two are in the group "ordinates"

        self.abscissae.fromkeys(self.infile['/distributions/rhoDist/abscissae'].keys(),0)
        ordinate = self.infile['/distributions/rhoDist/ordinate'].value

        #ADDING DIFFERENT ION SPECIES
        if ordinate.shape[-1]>25:
            diff=ordinate.shape[-1]-25
            n_ions_more=diff/2
            self.name_dict['ctor_el']+= n_ions_more
            self.name_dict['ctor_i1']+= n_ions_more              
            for el in range(n_ions_more):
                k='pi'+str(el+2)
                self.name_dict[k]=24+el
                k2='ctor_i'+str(el+2)
                self.name_dict[k2]=26+el+1

        
        for key in self.abscissae:
            self.abscissae[key]=self.infile[tree_path+'/abscissae/'+str(key)]
        self.rho = np.linspace(0,1,len(self.infile['/distributions/rhoDist/abscissae/dim1'])-1)


        #self. ordinate structure WILL BECOME:
        #(1,1,....,injector, time, rho, type of distribution)
        self.slices = ordinate.reshape(ordinate.shape[-4], ordinate.shape[-3], ordinate.shape[-2], ordinate.shape[-1])
        self.n_inj = self.slices.shape[0]
        self.dim_num = len(self.infile['/distributions/rhoDist/ordinates/'].keys())
        self.dim_num = int(self.dim_num/2.) #this because name and unit are double
        self.lab = np.array([], dtype='S32')
        self.uni = np.array([], dtype='S8')
        for i in range(self.dim_num):
            self.lab = np.append(self.lab, self.infile[tree_path+'ordinates/name_'+'{:06d}'.format(i+1)].value)
            self.uni = np.append(self.uni, self.infile[tree_path+'ordinates/unit_'+'{:06d}'.format(i+1)].value)

        #in infile['species/testParticle/origin'] there is an array with the ordered set of beams
        self._h5origins = self.infile['species/testParticle/origin'].value
        self._calc_originbeam()
            
    def TCV_plot_profs(self, *args):
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

    def plot_profs(self, inlist):
        """
        Method to plot quantities without different ion species (i.e. j,p,ecc.)
        """
        x = self.rho
        y_ind = np.array([(self.name_dict[t]-1) for t in inlist])
        y = np.zeros((self.n_inj, np.size(y_ind), np.size(x)))

        #(injector, time, rho, type of distribution)
        for t in range(np.size(y_ind)):
            y[:,t,:] = np.array([self.slices[:,0,:,y_ind[t]]])

        #y = np.array(self.slices[:,y_ind[:]])
        n_el = np.size(y_ind)
        ylabels = np.array(self.lab[y_ind[:]])
        yunits  = np.array(self.uni[y_ind[:]])
  
        n_row = 1
        n_col = 3
        if n_el < 3:
            n_col = n_el
        elif n_el > 3:
            n_row = 2
        elif n_el > 6:
            n_row = 3

        plt.figure()
        for i in range(n_el):
            plt.subplot(n_row,n_col,i+1)
            for j in range(self.n_inj):
                plt.plot(x, y[j,i,:], color=colours[0])
            plt.xlabel('rho')
            plt.ylabel(ylabels[i]+' '+yunits[i])

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
        if 1>2:
            self.beamlabel=['1','2','3','4','5','6','7','8','9','10',\
                            '13','14','NNB_U','NNB_L']
            self.beamlabel_num=['1','2','3','4','5','6','7','8','9','10',\
                                '13','14','99','101']
            self.beamorigin = np.array([45,46,47,48,133,134,135,136,\
                               221,222,223,224,309,310,311,312,\
                               3637,3638,3639,3640,5253,5254,5255,5256,\
                               3031,3032], dtype=int)
            self.beamlabel_no910 = ['1','2','3','4','5','6','7','8',\
                            '13','14','NNB_U','NNB_L']
        else:
            self.beamlabel=['1','2','3','4','5','6','9','10',\
                            '13','14','NNB_U','NNB_L']
            self.beamlabel_num=['1','2','3','4','5','6','9','10',\
                                '13','14','99','101']
            self.beamorigin = np.array([45,46,47,48,133,134,135,136,\
                               221,222,223,224,\
                               3637,3638,3639,3640,5253,5254,5255,5256,\
                               3031,3032], dtype=int)
            self.beamlabel_no910 = ['1','2','3','4','5','6','7','8',\
                            '13','14','NNB_U','NNB_L']            
        #orderedorigins is the array where you can retrieve the correct beam used in the distributions
        # so in orderedorigins we will find e.g. [13,7,9,11,15,4,8,9,...]
        # and you get data of beam 1 using orderedorigins[0], etc.
        self.orderedorigins=np.zeros(len(self._h5origins), dtype=int)
        j=0
        for i,el in enumerate(self.beamorigin):
            ind = np.where(self._h5origins==el)
            if not ind[0]:
                continue
            self.orderedorigins[j] = int(ind[0])
            j=j+1

    def calculate_scalar(self):
        """
        Method to calculate scalar quantities: total power to ions, to electrons,
        total current induced, total angular momentum
        """
        
        self._evaluate_shellArea()
        self._evaluate_shellVol ()
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
            self.calculate_scalar()
                

        for ii, el in enumerate(self.i_beams):
            print "Current from beam ", self.beamlabel[ii], " is ", el*1e-3, " kA"
        print ""
        for ii, el in enumerate(self.pe_beams):
            print "Pe from beam ", self.beamlabel[ii], " is ", el*1e-6, " MW"
        print ""
        for ii, el in enumerate(self.pi_beams):
            print "Pi from beam ", self.beamlabel[ii], " is ", el*1e-6, " MW"
        print ""
        for ii, el in enumerate(self.tor_beams):
            print "Tor from beam ", self.beamlabel[ii], " is ", el, " Nm"
            
        print " "
        print "Total current induced    ", self.I_tot*1e-3, " kA"
        print "Total power to electrons ", self.pe*1e-6, " MW"
        print "Total power to ions      ", self.pi1*1e-6, " MW"
        print "Total power delivered    ", (self.pi1+self.pe)*1e-6, " MW"
        print "Total torque ", self.tor, " Nm"

                
    def store_scal_to_ascii(self):
	"""
	Method to store data scalars (all grouped in a single file) to ascii
	"""

        try:
            np.mean(self.i_beams)
        except:
            self.calculate_scalar()
        num_col=4
        label_line='IND'+" "*7+'I (kA)'+" "*5+'Pe (MW)'+" "*5+'Pi (MW)'+" "*5
        data=np.zeros((len(self.beamlabel_num), num_col), dtype=float)
        data[:,0] = self.beamlabel_num
        data[:,1] = self.i_beams
        data[:,2] = self.pe_beams
        data[:,3] = self.pi_beams

        np.savetxt('scalar_data.dat', data, fmt='%.4e', header=label_line)

    def store_dis_to_ascii(self):
	"""
	Method to store data distribution (one for each beam) to ASCII
	"""
        header=''
        self._calc_originbeam  ()

        for oo,eloo in enumerate(self.lab):
            header+=self.lab[oo]+self.uni[oo]+'  '
        for ii, el in enumerate([1,2,3,4,5,6,7,8,9,10,13,14]):
            fname = str(el)+'.dat'
            data=self.slices[self.orderedorigins[2*ii],0,:,:]+self.slices[self.orderedorigins[2*ii+1],0,:,:]
            np.savetxt(fname, data, fmt='%.4e', header=header)

    def group_beams(self):
	"""
	Method to group data distribution for beam type (PPERP, PPAR, NNB)
	"""
        self._calc_originbeam  ()
                
        self.data_PPERP = np.zeros((len(self.rho), len(self.lab)),dtype=float)
        self.data_PPAR  = np.zeros((len(self.rho), len(self.lab)),dtype=float)
        self.data_NNB   = np.zeros((len(self.rho), len(self.lab)),dtype=float)

        for ii, el in enumerate(self.beamlabel):
        #for ii, el in enumerate(['NNB_L', 'NNB_U']):
            try:
                if el in [9,10]:
                    ind1=self.orderedorigins[2*ii]
                    ind2=self.orderedorigins[2*ii+1]
                    self.data_PPAR  = self.data_PPAR+self.slices[ind1,0,:,:]+self.slices[ind2,0,:,:]

                elif el=='NNB_L' or el=='NNB_U':
                    if el=='NNB_L':
                        tmp=-2
                    else:
                        tmp=-1
                    ind1=self.orderedorigins[tmp]
                    self.data_NNB   = self.data_NNB+ self.slices[ind1, 0,:,:]
                else:
                    ind1=self.orderedorigins[2*ii]
                    ind2=self.orderedorigins[2*ii+1]
                    self.data_PPERP = self.data_PPERP+self.slices[ind1,0,:,:]+self.slices[ind2,0,:,:]
            except:
                continue

    def group_beams_w78(self):
	"""
	Method to group data distribution for beam type (PPERP, PPAR, NNB)
	"""
        self._calc_originbeam  ()
                
        self.data_PPERP = np.zeros((len(self.rho), len(self.lab)),dtype=float)
        self.data_PPAR  = np.zeros((len(self.rho), len(self.lab)),dtype=float)
        self.data_NNB   = np.zeros((len(self.rho), len(self.lab)),dtype=float)

        for ii, el in enumerate(self.beamlabel_no910):
        #for ii, el in enumerate(['NNB_L', 'NNB_U']):
            try:
                if ii+1 in [7,8,9,10]:
                    ind1=self.orderedorigins[2*ii]
                    ind2=self.orderedorigins[2*ii+1]
                    self.data_PPAR  = self.data_PPAR+self.slices[ind1,0,:,:]+self.slices[ind2,0,:,:]

                elif el=='NNB_L' or el=='NNB_U':
                    if el=='NNB_L':
                        tmp=-2
                    else:
                        tmp=-1
                    ind1=self.orderedorigins[tmp]
                    self.data_NNB   = self.data_NNB+ self.slices[ind1, 0,:,:]
                else:
                    ind1=self.orderedorigins[2*ii]
                    ind2=self.orderedorigins[2*ii+1]
                    self.data_PPERP = self.data_PPERP+self.slices[ind1,0,:,:]+self.slices[ind2,0,:,:]
            except:
                continue
        
    def store_groupdis_to_ascii(self):
	"""
	Method to store data distribution (one for each BEAM TYPE) to ASCII
	"""
        try:
            self.data_PPAR.mean()
        except:
            self.group_beams()
        header=''

        for oo,eloo in enumerate(self.lab):
            header+=self.lab[oo]+self.uni[oo]+'  '
                
        np.savetxt('PPERP.dat', self.data_PPERP, fmt='%.4e', header=header)
        np.savetxt('PPAR.dat' , self.data_PPAR , fmt='%.4e', header=header)
        np.savetxt('NNB.dat'  , self.data_NNB  , fmt='%.4e', header=header)
        


    def plot_groupcurrent(self):
        """
        method to plot the data produced with group_beams 
        """
        try:
            self.data_PPAR.mean()
        except:
            self.group_beams()

        i_ppar = self.data_PPAR [:,2]*1e-3
        i_pper = self.data_PPERP[:,2]*1e-3
        i_nnb  = self.data_NNB  [:,2]*1e-3


        i_t_ppar = np.abs(self.data_PPAR [:,12])*1e-3
        i_t_pper = np.abs(self.data_PPERP[:,12])*1e-3
        i_t_nnb  = np.abs(self.data_NNB  [:,12])*1e-3
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

        thn_ppar=self.data_PPAR  [:,18]
        thn_pper=self.data_PPERP [:,18]
        thn_nnb=self.data_NNB    [:,18]

        press_ppar = self.data_PPAR [:,13]*1e-3
        press_pper = self.data_PPERP[:,13]*1e-3
        press_nnb  = self.data_NNB  [:,13]*1e-3
        p_tot      = press_ppar+press_pper+press_nnb

        tor_i_ppar = self.data_PPAR [:,5]+self.data_PPAR [:,24]
        tor_i_pper = self.data_PPERP[:,5]+self.data_PPERP[:,24]
        tor_i_nnb  = self.data_NNB  [:,5]+self.data_NNB  [:,24]
        tor_i_tot  = tor_i_ppar+tor_i_pper+tor_i_nnb
        
        labels=['P-T','P-P','N-NB']
        labels2=['P-T','P-P','N-NB','TOT']
#        plot_article(4,[self.rho, i_ppar, i_pper, i_nnb, i_tot],labels2,r'$\rho$', 'j (kA/$m^2$)', self.infile_n)
#        plot_article(4,[self.rho, n_ppar, n_pper, n_nnb, n_tot],labels2,r'$\rho$', 'n ($m^{-3}$)', self.infile_n)
#        plot_article(4,[self.rho, pe_ppar, pe_pper, pe_nnb, pe_tot],labels2,r'$\rho$', '$P_e$ (kW/$m^3$)', self.infile_n)
#        plot_article(4,[self.rho, pi_ppar, pi_pper, pi_nnb, pi_tot],labels2,r'$\rho$', '$P_i$ (kW/$m^3$)', self.infile_n)
#        plot_article(4,[self.rho, press_ppar, press_pper, press_nnb, p_tot],labels2,r'$\rho$', '$p$ (kPa)', self.infile_n)
#        plot_article(4,[self.rho, tor_i_ppar, tor_i_pper, tor_i_nnb, tor_i_tot],labels2,r'$\rho$', 'Torque density (N $m^{-2}$)', self.infile_n)

#        plot_article(3,[self.rho, tor_el_ppar, tor_el_pper, tor_el_nnb],labels,r'$\rho$', 'Torque to electrons (N $m^{-2}$)')
        plot_article(4,[self.rho, i_t_ppar, i_t_pper, i_t_nnb, i_tot],labels2,r'$\rho$', 'Shielded current (kA/$m^{2}$)', self.infile_n)

        #plot_article(3,[self.rho, thn_ppar, thn_pper, thn_nnb],labels,r'$\rho$', 'Th. n (1/$m^3$)')

        plt.show()



    def sum_all(self):
        newprofs = np.sum(self.slices, axis=0)
        newprofs = np.sum(newprofs, axis=0)
        self.Itot_prof = newprofs[:,12]
        self.ntot_prof = newprofs[:,0]
        self.petot_prof = newprofs[:,21]
        self.pitot_prof = newprofs[:,22]
        if self.slices.shape[-1]>25:
            self.pimptot_prof = newprofs[:,23]
        


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


def plot_article(n_lines, data, data_labels, xlabel, ylabel, suptitle):
        #=====================================================================================
        # SET TEXT FONT AND SIZE
        #=====================================================================================
        plt.rc('font', family='serif', serif='Palatino')
        #plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=====================================================================================
        col=['r','g','b','k']
        fig=plt.figure()
        fig.suptitle(suptitle)
        ax=fig.add_subplot(111)
        for i in range(n_lines):
            ax.plot(data[0], data[i+1], label=str(data_labels[i]), linewidth=3, color=col[i])
        ax.plot(data[0], np.zeros(len(data[0])), 'k',linewidth=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid('on')
        ax.set_ylim([0,250])
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
        ax.legend(loc='best')
