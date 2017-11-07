# ASCOT-Python
Contains the classes and methods built in python to create input and analyse output for ASCOT

ascot_Bfield.py - contains two classes:
    - Bfield_out: to import the equilibrium data computed by ascot, stored in the h5 file
    - Bfield_in : to take the input from an eqdsk file, and eventually build and write the input.magn_header and input.magn_bkg   needed for ASCOT and BBNBI
    
    
ascot_profs.py - contains the following classes:
    - profiles (SUPERCLASS): this class contains the methods shared all over all the classes (extrapolation over rho=1, plot, writing to files
    - h5_profiles: reads from h5 file (both BBNBI or ASCOT)
    - dat_profiles: reads from a list of ascii files in the form (rho, data), one for each quantity needed)
    - matlab_profiles: reads from a matlab file, usually produced with METIS code
    - ufiles: reads from a list of ufiles
    - SA_datfiles: reads from the data from CRONOS usually in the EU repository
   
   
ascot_particles.py - 


ascot_distributions.py - 
