"""
matteo.vallar@igi.cnr.it
Classes to read magnetic bkg input
"""
import numpy as np
      
def read_bkg(fname='input.magn_bkg'):
    """
    """
    f = open(fname, 'r'); d=f.readlines(); f.close()
    bkg = {}
    line = d[0].split()
    bkg['phi0'] = float(line[0])
    bkg['nsector'] = int(line[1])
    bkg['nphi_per_sector'] = int(line[2])
    bkg['ncoil'] = int(line[3])
    bkg['zero_at_coil'] = int(line[4])
    
    line = d[1].split()
    Rmin = float(line[0]); Rmax = float(line[1]); nR = int(line[2])
    bkg['R'] = np.linspace(Rmin, Rmax, nR)
    line = d[2].split()
    zmin = float(line[0]); zmax = float(line[1]); nz = int(line[2])
    bkg['z'] = np.linspace(zmin, zmax, nz)
    
    if bkg['nsector']==0:
        ind_psi = 5
    else:
        return()
        
    for rr in ['psi', 'BR', 'Bphi', 'Bz']:
        arr=np.array([])
        for i in range(int(np.floor(nR*nz/4.)+1)):
            for jj, el in enumerate(d[ind_psi+i].split()):
                arr = np.append(arr, float(el))
        bkg[rr] = np.reshape(arr, (nR, nz))
        ind_psi += i

    return bkg

def read_header(fname='input.magn_header'):
    """
    """
    f = open(fname, 'r'); d=f.readlines(); f.close()
    hdr = {}
    line = d[0].split()
    hdr['nSHOT'] = line[0]
    hdr['tSHOT'] = float(line[1])
    hdr['modflg'] = int(line[2])

    line = d[1].split()
    hdr['devnam'] = line[0]

    line = d[2].split()
    hdr['FPPkat'] = line[0]
    hdr['IpiFPP'] = line[1]

    line = d[3].split()
    nPF = int(line[0])
    hdr['PFxx'] = np.zeros(nPF)    
    hdr['RPFx'] = np.zeros(nPF)
    hdr['zPFx'] = np.zeros(nPF)
    line = d[4].split()
    for j in range(nPF):
        hdr['PFxx'][j] = line[j]

    line = d[5].split()
    for j in range(nPF):
        hdr['RPFx'][j] = line[j]

    line = d[6].split()
    for j in range(nPF):
        hdr['zPFx'][j] = line[j]
        
    hdr['SSQ'] = np.array([])
    line = d[7].split(); i=0;
    while np.size(line)>1:
        i=i+1    
        hdr['SSQ'] = np.append(hdr['SSQ'], (float(line[j]) for j in [0,1,2,3]))
        line = d[7+i].split()

    hdr['rhoPF'] = int(line[0])
    ind_arr = 7+i
    for rr in ('PFL','Vol','Area','Qpl'):
        arr=np.array([])
        for i in range(int(np.floor(hdr['rhoPF']/4.)+1)):
            for jj, el in enumerate(d[ind_arr+i].split()):
                arr = np.append(arr, float(el))
        hdr[rr] = arr
        ind_arr += i

    return hdr
