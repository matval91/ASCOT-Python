"""
Class for magnetic field
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib
#from matplotlib.figure import Figure

class Bfield:
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



    
    def __init__(self, infile):
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

        self.vardict = {}
        #store data with the correct labels
        
        for k in self.labdict.keys():
            self.vardict[k] = infile[self.labdict[k]]

        self.nrho = self.vardict['psi_rho'].shape[0]
        self.rho = np.zeros(self.nrho)
        self.rho = self.vardict['psi_rho'][:]
        self.nthe = self.vardict['theta'].shape[0]
        self.the = np.zeros(self.nthe)
        self.the = self.vardict['theta'][:]


        self.rmin = self.vardict["rmin"][:]
        self.zmin = self.vardict["zmin"][:]
        self.rmax = self.vardict["rmax"][:]
        self.zmax = self.vardict["zmax"][:]

        #these nR and nZ are from the group /bfield/, the ones for /boozer/ are on rho and theta (spacing of grid)
        self.nR = self.vardict['nr'][:]        
        self.nZ = self.vardict['nz'][:]
        self.R = np.linspace(self.rmin, self.rmax, self.nR)
        self.Z = np.linspace(self.zmin, self.zmax, self.nZ)
        
        # doing plot to check what we have
        #self.checkplot()


    def plotB(self):
        """
        Method to plot the magnetic field
        bphi, btheta are in 2D
        """
        
        bphi   = np.zeros((self.nR, self.nZ), dtype = float)
        bphi   = self.vardict["bphi"][:]
        btheta = np.zeros((self.nthe, self.nrho), dtype = float)
        btheta = self.vardict["btheta"][:]
        #btheta = np.transpose(btheta)
        print(np.shape(btheta))
        babs   = np.zeros((self.nrho, self.nthe), dtype = float)
        babs   = self.vardict["Babs"][:]
        
        n_row = 1
        n_col = 3
        plt.subplot(n_row,n_col,1)
        plt.title("B PHI")
        r,z = self.R, self.Z
        CS = plt.contourf(r,z,bphi, self.nrho)
        #CB = plt.colorbar(CS, shrink=0.8, extend='both')
        plt.axis([self.rmin, self.rmax, self.zmin, self.zmax])
        plt.xlabel("R")
        plt.ylabel("Z")
        
        plt.subplot(n_row,n_col,2, projection = 'polar')
        plt.title("B THETA")
        rho,theta = self.rho, self.the
        print rho.shape, theta.shape, np.shape(btheta)
        CS = plt.contourf(theta,rho, btheta, self.nrho)
        #CB = plt.colorbar(CS, shrink=0.8, extend='both')
        #plt.xlabel("R")
        #plt.ylabel("Z")
        
        plt.subplot(n_row,n_col,3, projection = 'polar')
        plt.title("B ABS")
        rho,theta = self.rho, self.the
        CS = plt.contourf(rho, theta, babs, self.nrho)
        #CB = plt.colorbar(CS, shrink=0.8, extend='both')
        #plt.xlabel("R")
        #plt.ylabel("Z")

        plt.show()	


    def checkplot(self):
        """
        Method to plot the values and check the magnetic field we are looking at
        
        """

        n_row = 1
        n_col = 3
        for e,val in enumerate(["psi_2D","q","g"]):

            plt.subplot(n_row,n_col,e+1)
            plt.title(val)

            #plot for mag surfaces 2D
            if val == "psi_2D":
                r,z = self.vardict['r'], self.vardict['z']
                CS = plt.contourf(r,z,self.vardict[val], self.nrho)
                CB = plt.colorbar(CS, shrink=0.8, extend='both')
                plt.xlabel("R")
                plt.ylabel("Z")
            
            else: #plot of q and g
                plt.plot(self.rho, self.vardict[val])
                plt.xlabel('rho')
                plt.ylabel(val)

        plt.show()
