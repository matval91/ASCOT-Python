"""
Class for profiles
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

colours = ['k','b','r','c','g']

class profiles:
    """
    Class for handling the profiles data
    
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



    
    def __init__(self, infile):
        #defining dictionary for handling data from h5file
        self.labdict = {"rho":"/plasma/1d/rho",\
                        "ne":"/plasma/1d/ne","ni":"/plasma/1d/ni",\
                        "te":"/plasma/1d/te","ti":"/plasma/1d/ti",\
                        "vtor":"/plasma/1d/vtor","zeff":"/plasma/1d/zeff",\
                        "a_ions":"/plasma/anum", "znum":"/plasma/znum"\
                        }

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
        Method to plot the values and check the profiles data we have
        
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


