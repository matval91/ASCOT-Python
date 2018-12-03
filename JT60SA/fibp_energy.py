"""
Script to join different FIBP profiles (w 1e6 particles for each set of beams)
"""
import ascot_distributions
from ascot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt
colours = ['k', 'r', 'b']
shot = 5
folder = '/home/vallar/ASCOT/runs/JT60SA/'
    
folder += '005/'
bfield_name = 'ascot_005065.h5'
dd = {
    '500':{
        'fname':'bbnbi_005050.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[]},
    '350':{
        'fname':'bbnbi_005058.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[]},
    '250':{
        'fname':'bbnbi_005056.h5',
        'dist_obj': None,
        'fibp': [],
        'rho':[]},
}

dist = ascot_distributions.SA_1d(folder+bfield_name)
rho = np.linspace(0, 1, num=len(dist._volumes), dtype=float)
for el in dd:
    bb_fname = folder+dd[el]['fname']
    dd[el]['rho']=rho
    dist.SA_calc_FIBP(bb_fname=bb_fname)        
    dd[el]['fibp'] = dist.fibp


data = np.array([rho, dd['500']['fibp'].T, dd['350']['fibp'].T, dd['250']['fibp'].T])

data_labels = ['500','350','250']
plot_article(3,data, data_labels=data_labels,xlabel=r'$\rho$', ylabel=r'Fast ion birth profile $1/(s\cdot m^3)$', ylim=[0, 0.3e19], col=colours)

