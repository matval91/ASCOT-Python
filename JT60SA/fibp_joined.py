"""
Script to join different FIBP profiles (w 1e6 particles for each set of beams)
"""
import a4py.classes.distributions as ascot_distributions
from utils.plot_utils import plot_article
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

shot = input('Scenario? ')
folder = '/home/vallar/ASCOT/runs/JT60SA/'
shot=int(shot)
if shot==2:
    folder += '002/'
    bfield_name = 'ascot_002034.h5'
    dd = {
        'NNB': {
            'fname':'bbnbi_002046.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'NNB'},
        'PT': {
            'fname':'bbnbi_002048.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PT'},
        'PP': {
            'fname':'bbnbi_002050.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PP'},
        'tot':{
            'fname':'bbnbi_002044.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[]}
    }
elif shot==3:
    folder += '003/'
    bfield_name = 'ascot_003001.h5'
    dd = {
        'NNB': {
            'fname':'bbnbi_003004.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'NNB'},
        'PT': {
            'fname':'bbnbi_003006.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PT'},
        'PP': {
            'fname':'bbnbi_003008.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PP'},
        'tot':{
            'fname':'bbnbi_003002.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[]}
    }
elif shot==4:
    folder += '004/'
    bfield_name = 'ascot_004008.h5'
    dd = {
        'NNB': {
            'fname':'bbnbi_004021.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'NNB'},
        'PT': {
            'fname':'bbnbi_004028.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PT'},
        'PP': {
            'fname':'bbnbi_004022.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PP'},
        'tot':{
            'fname':'bbnbi_004023.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[]}
    }
elif shot==5:
    folder += '005/'
    bfield_name = 'ascot_005065.h5'
    dd = {
        'NNB': {
            'fname':'bbnbi_005050.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'NNB'},
        'PP': {
            'fname':'bbnbi_005054.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PT'},
        'PT': {
            'fname':'bbnbi_005052.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[],
            'label': 'PP'},
        'tot':{
            'fname':'bbnbi_005048.h5',
            'dist_obj': None,
            'fibp': [],
            'rho':[]}
    }

dist = ascot_distributions.SA_1d(folder+bfield_name)
rho = np.linspace(0, 1, num=len(dist._volumes), dtype=float)
dd['tot']['fibp'] = np.zeros(len(rho))
for el in dd:
    bb_fname = folder+dd[el]['fname']
    dd[el]['rho']=rho
    if el != 'tot':
        dist.SA_calc_FIBP(bb_fname=bb_fname)        
        dd[el]['fibp'] = dist.fibp
        dd['tot']['fibp']+=dist.fibp
    fibp = dd[el]['fibp']


data = np.array([rho, dd['tot']['fibp'].T, dd['NNB']['fibp'].T, dd['PT']['fibp'].T, dd['PP']['fibp'].T])

data_labels = ['TOT','NNB', 'PT', 'PP']
plot_article(4,data, data_labels=data_labels,xlabel=r'$\rho$', ylabel=r'Fast ion birth profile $1/(s\cdot m^3)$',ylim=[0, 3.e19])

