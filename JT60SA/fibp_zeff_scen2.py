"""
Script to join different FIBP profiles (w 1e6 particles for each set of beams)
"""
import a4py.classes.distributions as ascot_distributions
from utils.plot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt

shot = 2
folder='/home/vallar/ASCOT/runs/SA_002/zeff_scan/fullpower'

bfield_name = '/ascot_dist.h5'
dd = {
    'Ref.':{
        'fname':'/bbnbi.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[], 'lc':'k'},
    'Core':{
        'fname':'_core/bbnbi.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[], 'lc':'b'},
    'Edge':{
        'fname':'_edge/bbnbi.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[], 'lc':'r'},
}

dist = ascot_distributions.SA_1d(folder+bfield_name)
rho = np.linspace(0, 1, num=len(dist._volumes), dtype=float)
col=[]

for el in dd:
    bb_fname = folder+dd[el]['fname']
    dd[el]['rho']=rho
    dist.SA_calc_FIBP(bb_fname=bb_fname)        
    dd[el]['fibp'] = dist.fibp


data = np.array([rho, dd['Ref.']['fibp'].T, dd['Core']['fibp'].T, dd['Edge']['fibp'].T])

data_labels = ['Ref.', 'Acc.', 'Edge']
col= [dd['Ref.']['lc'],dd['Core']['lc'],dd['Edge']['lc']]
plot_article(3,data, data_labels=data_labels,xlabel=r'$\rho$', ylabel=r'Fast ion birth profile $1/(s\cdot m^3)$', ylim=[0, 3.e19], col=col)

