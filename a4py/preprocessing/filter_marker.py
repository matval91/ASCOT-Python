
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
def filter_marker(input_fname='input.particles', \
                  fname_out='input.particles_filt',\
                  minrho=0.8, minxi=-1., maxxi=1., sign=1, max_markers=30000, maxrho=5.):

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
    | markers_towrite (array) : matrix with data of markers selected
    | markers         (array) : matrix with read data
    | pitch           (array) : pitch of data to be written
    """

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
        elif 'magnetic' in tmpl:
            if 'toroidal' in tmpl:
                indBphi = ind-ind_countrho-1
            elif 'vertical' in tmpl:
                indBz = ind-ind_countrho-1
            elif 'radial' in tmpl:
                indBr = ind-ind_countrho-1

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
    
    modB = np.sqrt(markers[:, indBphi]**2+markers[:, indBz]**2+markers[:, indBr]**2)
    vdotB = markers[:, indBphi]*markers[:, indvphi]+markers[:, indBz]*markers[:, indvz]+markers[:, indBr]*markers[:, indvr]
    pitch = vdotB/(modB*vtot)
    #vvparallel = vdotB/modB**2*[markers[:, indBphi], markers[:, indBz], markers[:, indBr]]
    #vparallel = np.sqrt(np.power(vvparallel[0,:],2)+np.power(vvparallel[1,:],2)+np.power(vvparallel[2,:],2))
    #pitch = sign*vparallel/vtot
    
    indnew = np.where(np.logical_and(np.logical_and(markers[:,indrho] > minrho, markers[:,indrho]<maxrho),\
                                     np.logical_and(pitch>minxi, pitch<maxxi)))[0]
    _markers_towrite = markers[indnew,:]
    if np.shape(indnew)[0]<max_markers: max_markers=np.size(indnew)
    max_markers=int(max_markers)
    markers_towrite = _markers_towrite[0:max_markers,:]
    n_newmarkers = max_markers

    tmp=header[ind_nmarkers].split("#")
    tmp[0] = str(n_newmarkers)
    header[ind_nmarkers] = " # ".join(tmp)
    header = "".join(header)
    fmt = ['%i','%7.6e','%i','%7.6e','%7.6e', '%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e','%7.6e']
    fmt[0] = '%i'; fmt[2] = '%i'; fmt[12] = '%i'; fmt[14] = '%i'
    np.savetxt(fname_out, markers_towrite, fmt=fmt,header=header, footer='#EOF', newline='\n', comments='')
    return markers, markers_towrite, pitch
