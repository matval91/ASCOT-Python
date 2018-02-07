"""
Class for distributions
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

colours = ['k','g','c','b','r']

class distribution:
    """
    Class for handling the profiles data
    
    DATA in h5 file (09/01/2017, ascot4)
    state can be both inistate or endstate

    /distributions/rhoDist   Group
    /distributions/rhoDist/abscissae Group
    /distributions/rhoDist/abscissae/dim1 Dataset {201}
    /distributions/rhoDist/abscissae/dim2 Dataset {2}
    /distributions/rhoDist/abscissae/dim3 Dataset {2}

    /distributions/rhoDist/ordinate Dataset {1, 1, 1, 1, 1, 200, 29}
    
    /distributions/rhoDist/ordinates Group
    /distributions/rhoDist/ordinates/name_000001 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000002 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000003 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000004 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000005 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000006 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000007 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000008 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000009 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000010 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000011 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000012 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000013 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000014 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000015 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000016 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000017 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000018 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000019 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000020 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000021 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000022 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000023 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000024 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000025 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000026 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000027 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000028 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/name_000029 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000001 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000002 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000003 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000004 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000005 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000006 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000007 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000008 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000009 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000010 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000011 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000012 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000013 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000014 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000015 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000016 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000017 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000018 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000019 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000020 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000021 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000022 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000023 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000024 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000025 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000026 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000027 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000028 Dataset {SCALAR}
    /distributions/rhoDist/ordinates/unit_000029 Dataset {SCALAR}
    
    /distributions/rhoDist/shellArea Dataset {200}
    /distributions/rhoDist/shellVolume Dataset {200}
    /distributions/rhoDist/variance Dataset {1, 1, 1, 1, 1, 200, 31}



    # name and unit for the /distributions/rhoDist/ordinates/
    # got with the command infile['/distributions/rhoDist/ordinates/'].keys()
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
    22 Power deposition to electrons W m^{-3}
    23 Power deposition to background species  1 W m^{-3}
    24 Power deposition to background species  2 W m^{-3}
    25 Power deposition to background species  3 W m^{-3}
    26 Collisional torque deposition to electrons N m^{-2}
    27 Collisional torque deposition to background species  1 N m^{-2}
    28 Collisional torque deposition to background species  2 N m^{-2}
    29 Collisional torque deposition to background species  3 N m^{-2}


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

    
    METHODS:
    __init__(self, infile) to store variables
    plot_CX (self)         to plot values relative to Charge-eXchange
    plot_profs(self,inlist) to plot values without different effect for different ions
        e.g. J, P, ...
    plot_profs_wions(self_inlist) same as plot_profs, but including ions
        e.g. p_dep, ...
    checkplot(self)        to plot values and check them

    """



    
    def __init__(self, infile):
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
                          'pel':22, 'pi1':23,'pi2':24,'pi3':25,\
                          'ctor_el':26, 'ctor_i1':27,'ctor_i2':28, 'ctor_i3':29\
                          }
        self.thlist   = ['n', 'e_den', 'th_n', 'th_e_n', 'th_torque' ]
        self.plist    = ['ptot', 'ppar', 'pperp']
        self.jlist    = ['jpar', 'jperp','jxB', 'tot_tor_j']
        self.pdeplist = ['pel', 'pi1', 'pi2', 'pi3']
        self.torlist  = ['ctor_el', 'ctor_i1', 'ctor_i2', 'ctor_i3'] 
        
        self.abscissae = {}
        # these two are in the group "ordinates"

        self.abscissae.fromkeys(infile['/distributions/rhoDist/abscissae'].keys(),0)
        self.ordinate = infile['/distributions/rhoDist/ordinate'].value
        #self.ordinates.fromkeys(infile['/distributions/rhoDist/ordinates'].keys(),0)
        
        for key in self.abscissae:
            self.abscissae[key]=infile[tree_path+'/abscissae/'+str(key)]
        self.rho = np.linspace(0,1,len(infile['/distributions/rhoDist/abscissae/dim1'])-1)
        # for key in self.ordinate:
        #    self.ordinate[key] =infile[tree_path+'/ordinate/ '+str(key)] 
        # for key in self.ordinates:
        #    self.ordinates[key]=infile[tree_path+'/ordinates/'+str(key)]
        self.slices = np.array(self.ordinate[0,0,0,0,0,:,:])

        self.dim_num = len(infile['/distributions/rhoDist/ordinates/'].keys())
        self.dim_num = int(self.dim_num/2.) #this because name and unit are double
        self.lab = np.array([], dtype='S32')
        self.uni = np.array([], dtype='S8')
        for i in range(self.dim_num):
            self.lab = np.append(self.lab, infile[tree_path+'ordinates/name_'+'{:06d}'.format(i+1)].value)
            self.uni = np.append(self.uni, infile[tree_path+'ordinates/unit_'+'{:06d}'.format(i+1)].value)
            

    def plot_profs(self, inlist):
        """
        Method to plot quantities without different ion species (i.e. j,p,ecc.)
        """
        x = self.rho
        y_ind = np.array([(self.name_dict[t]-1) for t in inlist])
        y = np.zeros((np.size(y_ind), np.size(x)))
        y = np.array([self.slices[:,y_ind[t]] for t in range(np.size(y_ind))])
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

        fig = plt.figure()
        for i in range(n_el):
            plt.subplot(n_row,n_col,i+1)
            plt.plot  (x, y[i,:], color=colours[0])
            plt.xlabel('rho')
            plt.ylabel(ylabels[i]+' '+yunits[i])

        plt.show()


    def plot_profs_wions(self,inlist):
        """
        Method to plot quantities with different ion species (i.e. CX,pow dep, ecc)
        suppose the list is in this format : [el, i1,i2,i3,...]
        """        
        x = self.rho
        y_ind = np.array([(self.name_dict[t]-1) for t in inlist])
        y = np.zeros((np.size(y_ind), np.size(x)))
        y = np.array([self.slices[:,y_ind[t]] for t in range(np.size(y_ind))])
        #y = np.array(self.slices[:,y_ind[:]])
        n_species = 3 #need to obtain this value automatically
        n_el = np.size(y_ind)-(n_species-1) # need to subtract the ohter ion species for the plot
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

        fig = plt.figure()
        for i in range(n_el):
            plt.subplot(n_row,n_col,i+1)
            plt.plot  (x, y[i,:], color=colours[0])
            plt.xlabel('rho')
            plt.ylabel(ylabels[i]+' '+yunits[i])
        for j in range(n_species):
            plt.plot(x,y[i+j,:],label='ion '+str(j+1), color=colours[j])

        plt.legend()
        plt.show()

        
    def plot_CX(self):
        """
        Method to plot values relative to charge-exchange
        """
        x  = self.rho
        y_ind = np.array([self.name_dict['CX_ionsource' ], self.name_dict['CX_ionensource' ],\
                 self.name_dict['CX_neutsource'], self.name_dict['CX_neutensource']])
        slices = np.array(self.ordinate[0,0,0,0,0,:,:])
        #y0,y1,y2,y3 = self.ordinate[0,0,0,0,0,:,y_ind[0]], self.ordinate[0,0,0,0,0,:,y_ind[1]],\
        #              self.ordinate[0,0,0,0,0,:,y_ind[2]], self.ordinate[0,0,0,0,0,:,y_ind[3]]
        y0, y1, y2, y3 = slices[:, y_ind[0]], slices[:, y_ind[1]],\
                         slices[:, y_ind[2]], slices[:, y_ind[3]]
        ylabels = np.array(self.lab[y_ind[t]-1] for t in y_ind)
        yunits  = np.array(self.uni[y_ind[:]])

        
        n_row, n_col = 1,2
        fig = plt.figure()
        fig.suptitle('Charge Exchange')

        plt.subplot(n_row,n_col,1)
        plt.plot(x, y0, label='Ion'     , color=colours[0])
        plt.plot(x, y1, label='Neutrals', color=colours[1])
        plt.xlabel('rho')
        plt.ylabel(ylabels[0]+' '+yunits[0])

        plt.subplot(n_row,n_col,2)
        plt.plot(x, y2, label='Ion'     , color=colours[0])
        plt.plot(x, y3, label='Neutrals', color=colours[1])
        plt.xlabel('rho')
        plt.ylabel(ylabels[2]+' '+yunits[2])
        #plt.ylabel('Energy source (J m^{-3})')
        plt.legend()

        plt.show()




    def checkplot(self):
        """
        Method to plot the values and check the magnetic field we are looking at
        
        """

        n_row = 3
        n_col = 2
        for e,val in enumerate(['ne','te','ni','ti','vtor','zeff']):

            plt.subplot(n_row,n_col,e+1)
            #plt.title(val)

            #plot for multiple species
            if val == 'ni':
                ylabel = val
                for i in np.arange(self.nspec):
                    plt.plot(self.rho, self.vardict[val][:,i], colours[i],\
                             label="A="+str(self.vardict['a_ions'][i]))

                plt.legend()
                
            elif val == 'te' or val == 'ti':
                y=self.vardict[val][:]
                y=y*0.001
                ylabel = str(val)+' [keV]'
            else: #plot ne,te,ti
                y = self.vardict[val][:]
                ylabel = str(val)

            plt.plot(self.rho,y)
            plt.xlabel('rho')
            plt.ylabel(ylabel)

        plt.show()
