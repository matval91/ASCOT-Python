"""
Script to join different profiles (w 1e6 particles for each set of beams)
"""
import ascot_distributions, ascot_particles
from ascot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

def interp1d(x,y,xnew):
    p = interp.interp1d(x,y)
    ynew = p(xnew)
    return ynew

shot=5
#shot = input('Scenario? ')
folder = '/home/vallar/ASCOT/runs/JT60SA/'
    
if shot==2:
    folder += '002/'
    dd = {
        '500': {
            'fname':'ascot_002036.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':500e3},
        '350': {
            'fname':'ascot_002044.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':350e3},
        '250': {
            'fname':'ascot_002042.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':250e3}
    }
elif shot==5:
    folder += '005/'
    dd = {
        'ref': {
            'fname':'ascot_005065.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pel':[],
            'pr':[], 'n':[],
            'pi':[], 'pi2':[],   
            'label': 'ref', 'lc':'k'},
        'acc': {
            'fname':'ascot_005077.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[],   'pi2':[],
            'label': 'acc', 'lc':'b'},
        'edge': {
            'fname':'ascot_005079.h5',
            'dist_obj': None,
            'rho':[],
            'j': [],
            'pr':[], 'n':[],
            'pel':[],
            'pi':[], 'pi2':[],
            'label': 'edge', 'lc':'r'}
    }
col=[]
for el in dd:
    fname = dd[el]['fname']
    p=ascot_particles.SA_iniend(folder+fname)
    p.endcondition()

    dist = ascot_distributions.SA_1d(folder+fname)
    dist._group_beams()
    dd[el]['rho'] = dist.rho
    dd[el]['j'] = dist.slices_summed[0,:,12]
    dd[el]['pel'] = dist.slices_summed[0,:,dist.peind-1]
    dd[el]['pi'] = dist.slices_summed[0,:,dist.peind]
    dd[el]['pr'] = dist.slices_summed[0,:,13]
    dd[el]['n'] = dist.slices_summed[0,:,0]
    dd[el]['pi2'] = dist.slices_summed[0,:,dist.peind+1]
    dist.calc_eff();
    dd[el]['eff'] = dist.eff
    #_,_,_,dd[el]['taus'] = dist._ecrit(E0=dd[el]['E'])
    
    #dist.print_scalars()
    dd[el]['gi'] = 1-dist.pe/(dist.pe+dist.pi1+dist.pi2)
    dist.print_scalars()


rho = dd['ref']['rho']


data_labels = ['Ref.', 'Acc.','Edge']
ylab = [r'j (kA/$m^2$)', r'$P_{el}\,(kW/m^3)$', r'$P_{D}\,(kW/m^3)$', r'p (kN/$m^2$)', r'n (1/$m^3$)', r'$P_{C}\,(kW/m^3)$']
factor = [-1e-3, 1e-3, 1e-3, 1e-3, 1, 1e-3]
col= [dd['ref']['lc'],dd['acc']['lc'],dd['edge']['lc']]
for i, el in enumerate(['j', 'pel', 'pi', 'pr','n', 'pi2']):
    #tot = factor[i]*dd['500'][el].T+factor[i]*dd['350'][el].T+factor[i]*dd['250'][el].T
    data = np.array([rho, factor[i]*dd['ref'][el].T, factor[i]*dd['acc'][el].T, factor[i]*dd['edge'][el].T])
    if el=='pel' or el=='pi' or el == 'pi2':
        if shot==4 or shot==5:
            ylim=[0,350]
        else:
            ylim=[0,200]
    else:
        ylim=0
    plot_article(3,data, data_labels=data_labels,xlabel=r'$\mathbf{\rho}$', ylabel=ylab[i], ylim=ylim, col=col)

#f=plt.figure(); ax=f.add_subplot(211); ax2=f.add_subplot(212)
#E=[250,350,500]
#eff = [dd['250']['eff'],dd['350']['eff'],dd['500']['eff']]
#gi = [dd['250']['gi'],dd['350']['gi'],dd['500']['gi']]
#taus = [dd['250']['taus'],dd['350']['taus'],dd['500']['taus']]

#ax.scatter(E, eff)
#ax.scatter(E, gi, marker='x')
#ax2.scatter(E, taus)
plt.show()
