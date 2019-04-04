import a4py.classes.orbits as orbits

from utils.plot_utils import _plot_1d, plot_article, limit_labels
import numpy as np
import matplotlib.pyplot as plt
import pickle

shot=2
#shot = input('Scenario? ')
folder = '/home/vallar/ASCOT/runs/JT60SA/'
if shot==2:

    folder = '/home/vallar/JT60-SA/002_lowden/impurity_effect/'
    fname = folder+'scen2_3e4.pkl'
elif shot==15:
    folder = '/home/vallar/ASCOT/runs/SA_005/reducedpower/reducedpower_allbeams'
    shots = ['Ref.','Core','Edge']
    colors = ['k', 'b', 'r']
    dd = {key:{} for key in shots}
    subkeys = ['fname', 'dist_obj', 'rho', 'j', 'pel', 'pr', 'n', 'pi', 'pi2', 'label', 'lc',\
               'dE']
    for ii, kk in enumerate(shots):
        dd[kk] = {key:[] for key in subkeys}
        dd[kk]['lc'] = colors[ii]
        dd[kk]['label']=kk
        dd[kk]['fname']=''
    dd['Ref.']['fname'] = '/orbits_ref.dat'
    dd['Core']['fname'] = '_acc/orbits_acc.dat'
    dd['Edge']['fname'] = '_edge/orbits_edge.dat'
    dd['Ref.']['fname_surf'] = '/ascot_dist.h5'
    dd['Core']['fname_surf'] = '_acc/ascot_dist.h5'
    dd['Edge']['fname_surf'] = '_edge/ascot_dist.h5'
    
col=[]
dict_data = pickle.load(open(fname, 'rb'))
colors = ['k', 'b', 'r']
dict_data['Ref.'][0]['lc']='k'; dict_data['Core'][0]['lc']='b'; dict_data['Edge'][0]['lc']='r'

if False:
    figr=plt.figure(); axrho=figr.add_subplot(111)
    figp=plt.figure(); axpitch=figp.add_subplot(111)
    figw=plt.figure(); axw=figw.add_subplot(211); axw_n = figw.add_subplot(212)
    figrp=plt.figure(); axrp=figrp.add_subplot(111)
    lim_banana_nnb=7.5*0.01
    for el in dict_data:
        rho = dict_data[el][0]['rho']
        col = dict_data[el][0]['lc']
        axrho.hist(rho, histtype='step', color=col, label=el, lw=2.3, density=True, bins=40)
        pitch = dict_data[el][0]['pitch']
        axpitch.hist(pitch, histtype='step', color=col,  label=el, lw=2.3, density=True, bins=40)
        axrp.scatter(rho, pitch, color=col, marker='o', alpha=0.3)
        w_banana = dict_data[el][0]['wbanana']
        indp = np.where(w_banana<lim_banana_nnb)[0]
        axw.hist(w_banana[indp]*100., histtype='step', label=el, color=col, lw=2.3, density=True, bins=40)
        indn = np.where(w_banana>lim_banana_nnb)[0]
        axw_n.hist(w_banana[indn]*100., histtype='step', label=el, color=col, lw=2.3, density=True, bins=40)

    limit_labels(axw, r'$w_{banana}$ [cm]'); limit_labels(axw_n, r'$w_{banana}$ [cm]')  
    limit_labels(axrho, r'$\rho$')  
    limit_labels(axpitch, r'$\xi$')      
    limit_labels(axrp, r'$\rho$',  r'$\xi$')  
    figr.tight_layout(); figw.tight_layout()
    figp.tight_layout()
plt.show()
