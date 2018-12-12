""" Filter markers

Script to filter the input.particles file (useful for wall-losses studies) 
basing on rho and xi

Parameters:
    | input_fname (str)   :  name of file to read (default input.particles)
    | fname_out   (str)   : name of file where to write (default input.particles_filt)
    | minrho      (float) : minimum rho (particles will be chosen after this rho value)\
                            (default is 0.8)
    | minxi       (float) : minimum pitch allowable (default is -1)
    | maxxi       (float) : maximum pitch allowable (default is 1)   
Arguments:
    | filter_marker: matrix with data of markers selected
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
def filter_marker(input_fname='input.particles', \
                  fname_out='input.particles_filt',\
                  minrho=0.8, minxi=-1., maxxi=1., sign=1):

    fin = open(input_fname,"r")
    lines=fin.readlines()
    for ind, ll in enumerate(lines):
        tmpl = ll.split()
        if 'fields' in tmpl:
            nfields = int(tmpl[0])
            ind_countrho = ind
        elif 'particles' in tmpl:
            nmarkers = int(tmpl[0])
            ind_nmarkers = ind
        elif 'flux' in tmpl:
            indrho = ind-ind_countrho-1
        elif 'velocity' in tmpl:
            if 'toroidal' in tmpl:
                indvphi = ind-ind_countrho-1
            elif 'vertical' in tmpl:
                indvz = ind-ind_countrho-1
            elif 'radial' in tmpl:
                indvr = ind-ind_countrho-1
        try:
            float(tmpl[1])
        except:
            continue
        ind_markerstart = ind
        break
    
    header = lines[0:ind_markerstart-1]
    markers = np.zeros((nmarkers, nfields))
    for ind, ll in enumerate(lines[ind_markerstart:-1]):
        tmp = ll.split()
        markers[ind,:] = tmp[:]
    vtot = np.sqrt(markers[:, indvphi]**2+markers[:, indvz]**2+markers[:, indvr]**2)
    pitch = sign*markers[:, indvphi]/vtot
    indnew = np.where(np.logical_and(markers[:,indrho] > minrho,\
                                     np.logical_and(pitch>minxi, pitch<maxxi)))[0]
    
    markers_towrite = markers[indnew,:]
    n_newmarkers = len(indnew)

    tmp=header[ind_nmarkers].split("#")
    tmp[0] = str(n_newmarkers)
    header[ind_nmarkers] = " # ".join(tmp)
    header = "".join(header)
    fmt = ['%i','%7.6e','%i','%7.6e','%7.6e', '%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e']
    fmt[0] = '%i'; fmt[2] = '%i'; fmt[12] = '%i'; fmt[14] = '%i'

    np.savetxt(fname_out, markers_towrite, fmt=fmt,header=header, footer='#EOF', newline='\n', comments='')
    return markers[indnew, indrho], pitch[indnew] 