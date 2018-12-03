"""
Script to join different FIBP profiles (w 1e6 particles for each set of beams)
"""
import ascot_distributions
from ascot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt

shot = 5
folder = '/home/vallar/ASCOT/runs/JT60SA/'
    
folder += '005/'
bfield_name = 'ascot_005065.h5'
dd = {
    'ref':{
        'fname':'bbnbi_005048.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[], 'lc':'k'},
    'acc':{
        'fname':'bbnbi_005060.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[], 'lc':'b'},
    'edge':{
        'fname':'bbnbi_005062.h5',
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


data = np.array([rho, dd['ref']['fibp'].T, dd['acc']['fibp'].T, dd['edge']['fibp'].T])

data_labels = ['Ref.', 'Acc.', 'Edge']
col= [dd['ref']['lc'],dd['acc']['lc'],dd['edge']['lc']]
plot_article(3,data, data_labels=data_labels,xlabel=r'$\rho$', ylabel=r'Fast ion birth profile $1/(s\cdot m^3)$', ylim=[0, 3.e19], col=col)

