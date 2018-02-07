"""
Class for profiles
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

colours = ['k','b','r','c','g']

class state:
    """
    Class for handling the profiles data
    
    DATA in h5 file (09/01/2017, ascot4)
    state can be both inistate or endstate
    
    /state/R              Dataset {10/Inf}
    /state/Rprt           Dataset {10/Inf}
    /state/cpuTime        Dataset {10/Inf}
    /state/dH             Dataset {10/Inf}
    /state/endcond        Dataset {10/Inf}
    /state/energy         Dataset {10/Inf}
    /state/fPoloidal      Dataset {10/Inf}
    /state/fToroidal      Dataset {10/Inf}
    /state/id             Dataset {10/Inf}
    /state/mhdPhase       Dataset {10/Inf}
    /state/npassing       Dataset {10/Inf}
    /state/ntrapped       Dataset {10/Inf}
    /state/origin         Dataset {10/Inf}
    /state/phi            Dataset {10/Inf}
    /state/phiprt         Dataset {10/Inf}
    /state/pitch          Dataset {10/Inf}
    /state/rho            Dataset {10/Inf}
    /state/species        Dataset {10/Inf}
    /state/time           Dataset {10/Inf}
    /state/vR             Dataset {10/Inf}
    /state/vphi           Dataset {10/Inf}
    /state/vz             Dataset {10/Inf}
    /state/wallTile       Dataset {10/Inf}
    /state/weight         Dataset {10/Inf}
    /state/z              Dataset {10/Inf}
    /state/zprt           Dataset {10/Inf}

    
    METHODS:
    __init__(self, infile) to store variables
    checkplot(self) to plot values and check them

    """



    
    def __init__(self, infile, statedef):
        if statedef != 'ini' && statedef != 'end':
            raise Error("Unknown state: "+str(statedef))
        
        #defining dictionary for handling data from h5file
        self.labdict = {}
        labels = ("R","Rprt","cputime","dH","endcond","energy","fPoloidal","fToroidal",\
                  "id","mhdPhase","npassing","ntrapped","origin","phi","phiprt",\
                  "pitch","rho","species","time","vR","vphi","vz",\
                  "wallTile","weight","z","zprt"\
                    
                    )
        for key in labels:
            label ='/'+state+'/'+ key
            self.labdict[key] = label
            
        self.vardict = {}
        #store data with the correct labels
        for k in self.labdict.keys():
            self.vardict[k] = infile[self.labdict[k]]

        self.rho = self.vardict['rho']
        self.nrho = self.rho.shape[0]
        self.nspec = len(self.vardict['a_ions'])
        if len(self.vardict['a_ions'])!=len(self.vardict['znum']):
            print "ERROR! array of A and Z don't have the same length"
            
        # plot to check what we have
        self.checkplot()




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


