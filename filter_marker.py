"""
Script to filter the input.particles file (useful for wall-losses studies)
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

def filter_marker(input_fname='input.particles', fname_out='input.particles_filt',minrho=0.8)

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

    minrho = 0.8
    indnew = np.where(markers[:,indrho] > minrho)[0]
    markers_towrite = markers[indnew,:]
    n_newmarkers = len(indnew)

    tmp=header[ind_nmarkers].split("#")
    tmp[0] = str(n_newmarkers)
    header[ind_nmarkers] = " # ".join(tmp)
    header = "".join(header)
    fmt = ['%i','%7.6e','%i','%7.6e','%7.6e', '%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e']
    fmt[0] = '%i'; fmt[2] = '%i'; fmt[12] = '%i'; fmt[14] = '%i'

    np.savetxt(fname_out, markers_towrite, fmt=fmt,header=header, footer='EOF#', newline='\n', comments='')