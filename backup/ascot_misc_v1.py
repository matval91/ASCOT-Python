"""
Class for misc
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt


class misc:
    """
    Class for handling the following misc specifications and plots
    Groups in h5 file (09/01/2017, ascot4)
    /wall                    Group
    /wall/2d                 Group
    /wall/2d/R               Dataset {116}
    /wall/2d/divFlag         Dataset {116}
    /wall/2d/segLen          Dataset {115}
    /wall/2d/z               Dataset {116}

    
    METHODS:
    __init__(self, infile) to store variables
    checkplot(self) to plot magnetic values and check them

    """



    
    def __init__(self, infile):
        #defining dictionary for handling data from h5file
        self.labdict = {"R_wall":"/wall/2d/R", "Z_wall":"/wall/2d/z",\
                        "divflag":"/wall/2d/divFlag", "segLen":"/wall/2d/segLen"\
                        }

        self.vardict = {}
        #store data with the correct labels
        
        for k in self.labdict.keys():
            self.vardict[k] = infile[self.labdict[k]]
        
        # doing plot to check what we have
        #self.checkplot()
	


    def checkplot(self):
        """
        Method to plot the values and check the magnetic field we are looking at
        
        """


        plt.title("2D wall")
        
        # plot for 2D wall
        plt.plot(self.vardict["R_wall"], self.vardict["Z_wall"])
        plt.xlabel("R")
        plt.ylabel("Z")

        plt.show()
