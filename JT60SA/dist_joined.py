"""
Script to join different profiles (w 1e6 particles for each set of beams)
"""
import ascot_distributions
from ascot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

def interp1d(x,y,xnew):
    p = interp.interp1d(x,y)
    ynew = p(xnew)
    return ynew




shot = input('Scenario? ')
folder = '/home/vallar/ASCOT/runs/JT60SA/'
    
if shot==2:
    folder += '002/'
    dd = {
        'NNB': {
            'fname':'ascot_002036.h5',
            'dist_obj': None,
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[],
            'label': 'NNB'},
        'PT': {
            'fname':'ascot_002038.h5',
            'dist_obj': None,
            'j': [],
            'pel':[],
            'pi':[],
            'pr':[],   'n':[],          
            'rho':[],
            'label': 'PT'},
        'PP': {
            'fname':'ascot_002040.h5',
            'dist_obj': None,
            'j': [],
            'pel':[],
            'pi':[], 
            'pr':[], 'n':[],            
            'rho':[],
            'label': 'PP'},
        'tot':{
            'fname':'ascot_002034.h5',
            'dist_obj': None,
            'fibp': [],
            'j': [],
            'pel':[],
            'pr':[],  'n':[],           
            'pi':[],   
            'rho':[]}
    }
elif shot==3:
    folder += '003/'
    dd = {
        'NNB': {
            'fname':'ascot_003003.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[],   
            'label': 'NNB'},
        'PT': {
            'fname':'ascot_003005.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   
            'label': 'PT'},
        'PP': {
            'fname':'ascot_003007.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   
            'label': 'PP'},
        'tot':{
            'fname':'ascot_003001.h5',
            'dist_obj': None,
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[],   
            'rho':[]}
    }
elif shot==4:
    folder += '004/'
    dd = {
        'NNB': {
            'fname':'ascot_004012.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[],   
            'label': 'NNB'},
        'PT': {
            'fname':'ascot_004019.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   
            'label': 'PT'},
        'PP': {
            'fname':'ascot_004014.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   
            'label': 'PP'},
        'tot':{
            'fname':'ascot_004018.h5',
            'dist_obj': None,
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[],   
            'rho':[]}
    }
elif shot==5:
    folder += '005/'
    dd = {
        'NNB': {
            'fname':'ascot_005067.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[],   
            'label': 'NNB'},
        'PT': {
            'fname':'ascot_005069.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   
            'label': 'PT'},
        'PP': {
            'fname':'ascot_005071.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   
            'label': 'PP'},
        'tot':{
            'fname':'ascot_005065.h5',
            'dist_obj': None,
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[],   
            'rho':[]}
    }

for el in dd:
    fname = dd[el]['fname']
    dist = ascot_distributions.SA_1d(folder+fname)
    dist.group_beams()
    dd[el]['rho'] = dist.rho
    if el=='NNB':
        dd[el]['j'] = dist.data_NNB[:,12]
        dd[el]['pel'] = dist.data_NNB[:,dist.peind-1]
        dd[el]['pi'] = dist.data_NNB[:,dist.peind]
        dd[el]['pr'] = dist.data_NNB[:,13]
        dd[el]['n'] = dist.data_NNB[:,0]
    elif el=='PT':
        dd[el]['j'] = dist.data_PPAR[:,12]
        dd[el]['pel'] = dist.data_PPAR[:,dist.peind-1]
        dd[el]['pi'] = dist.data_PPAR[:,dist.peind]
        dd[el]['pr'] = dist.data_PPAR[:,13]
        dd[el]['n'] = dist.data_PPAR[:,0]
    else:
        dd[el]['j'] = dist.data_PPERP[:,12]
        dd[el]['pel'] = dist.data_PPERP[:,dist.peind-1]
        dd[el]['pi'] = dist.data_PPERP[:,dist.peind]
        dd[el]['pr'] = dist.data_PPERP[:,13]        
        dd[el]['n'] = dist.data_PPERP[:,0]

rho = dd['NNB']['rho']


data_labels = ['TOT','PP', 'PT', 'NNB']
ylab = [r'j (kA/$m^2$)', r'$P_{el}\,(kW/m^3)$', r'$P_{i}\,(kW/m^3)$', r'p (kN/$m^2$)', r'n (1/$m^3$)']
factor = [-1e-3, 1e-3, 1e-3, 1e-3, 1]
for i, el in enumerate(['j', 'pel', 'pi', 'pr','n']):
    try:
        tot = factor[i]*dd['PP'][el].T+factor[i]*dd['PT'][el].T+factor[i]*dd['NNB'][el].T
        data = np.array([rho, tot, factor[i]*dd['PP'][el].T, factor[i]*dd['PT'][el].T, factor[i]*dd['NNB'][el].T])
    except:
        rho_old = dd['PT']['rho']
        dd['PT'][el] = interp1d(rho_old, dd['PT'][el],rho)
        tot = factor[i]*dd['PP'][el].T+factor[i]*dd['PT'][el].T+factor[i]*dd['NNB'][el].T
        data = np.array([rho, tot, factor[i]*dd['PP'][el].T, factor[i]*dd['PT'][el].T, factor[i]*dd['NNB'][el].T])

 
    if el=='pel' or el=='pi':
        if shot==4 or shot==5:
            ylim=[0,300]
        else:
            ylim=[0,200]
    else:
        ylim=0
    plot_article(4,data, data_labels=data_labels,xlabel=r'$\mathbf{\rho}$', ylabel=ylab[i], ylim=ylim)

