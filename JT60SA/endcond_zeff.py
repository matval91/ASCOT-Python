"""
Script to join different profiles (w 1e6 particles for each set of beams)
"""
import a4py.classes.distributions as ascot_distributions
import a4py.classes.particles as ascot_particles
import a4py.classes.distributions_2d as d2d
import a4py.classes.orbits as orbits

from utils.plot_utils import _plot_1d, plot_article, common_style
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

common_style()

shot=2
#shot = input('Scenario? ')
folder = '/home/vallar/ASCOT/runs/JT60SA/'
    
if shot==2:
    folder = '/home/vallar/ASCOT/runs/SA_002/zeff_scan/fullpower'
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
    dd['Ref.']['fname'] = '/ascot_dist.h5'
    dd['Core']['fname'] = '_core/ascot_dist.h5'
    dd['Edge']['fname'] = '_edge/ascot_dist.h5'
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
    dd['Ref.']['fname'] = '/ascot_dist.h5'
    dd['Core']['fname'] = '_acc/ascot_dist.h5'
    dd['Edge']['fname'] = '_edge/ascot_dist.h5'
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
    ind = np.where(p.data_e['endcond']==3)[0]
    dd[el]['particles'] = p
    dd[el]['dE'] = p.data_e['energy'][ind]
    dd[el]['phi'] = p.data_e['phi'][ind]
    dd[el]['rho'] = p.data_e['rho'][ind] 
    dd[el]['rhoini'] = p.data_i['rho'][ind]   
   
    theta = np.arctan2(p.data_e['z']-dd[el]['particles'].z0,p.data_e['R']-dd[el]['particles'].R0)
    dd[el]['theta'] = theta
    dd[el]['R'] = p.data_i['R'][ind]  
    p.endcondition()

f=plt.figure(); axrz=f.add_subplot(111)
for ss in shots:
    #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])
    data=dd[ss]['particles']
    x=np.cos(data.data_i['phi'])*data.data_i['R']
    y=np.sin(data.data_i['phi'])*data.data_i['R']
    axrz.scatter(x,y, marker='x', color=dd[ss]['lc'])
    #data=dd[ss]['theta']
    #axtheta.hist(data, histtype='step', color=dd[ss]['lc'], bins=20, lw=2.3, label=dd[ss]['label'])
plt.axis('equal')
f=plt.figure(); axr=f.add_subplot(111)
for ss in shots:
    #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])
    data=dd[ss]['rhoini']
    axr.hist(data, histtype='step', color=dd[ss]['lc'], bins=30, lw=2.3, label=dd[ss]['label'])
    #data=dd[ss]['theta']
    #axtheta.hist(data, histtype='step', color=dd[ss]['lc'], bins=20, lw=2.3, label=dd[ss]['label'])
axr.legend(loc='upper left')    
axr.set_xlabel(r'$\rho$')
axr.set_ylabel(r'Lost markers initial position (AU)', fontsize='x-large')
axr.grid('on')
axr.yaxis.set_ticks([])
f.tight_layout()

if True:
    f=plt.figure(); axe=f.add_subplot(111)
    bins = np.linspace(10, 90, 20)
    for ss in shots:
        #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])
        data=dd[ss]['dE']*1e-3
        axe.hist(data, histtype='step', color=dd[ss]['lc'], bins=30, lw=2.5, label=dd[ss]['label'])
    axe.legend(loc='upper center')
    axe.set_xlabel(r'E$_{final}$ [kev]')
    axe.set_ylabel(r'Lost markers (AU)')
    axe.grid('on')
    axe.yaxis.set_ticks([])
    f.tight_layout()

    f=plt.figure(); axphi=f.add_subplot(111)
    for ss in shots:
        #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])
        data=dd[ss]['phi']
        axphi.hist(data, histtype='step', color=dd[ss]['lc'], bins=30, lw=2.3, label=dd[ss]['label'])
        #data=dd[ss]['theta']
        #axtheta.hist(data, histtype='step', color=dd[ss]['lc'], bins=20, lw=2.3, label=dd[ss]['label'])

    axphi.legend(loc='lower center')    

    axphi.set_xlabel(r'$\phi$')
    axphi.set_ylabel(r'Lost markers (AU)')
    axphi.grid('on')
    axphi.yaxis.set_ticks([])
    axphi.xaxis.set_ticks([-np.pi, -0.5*np.pi, 0., 0.5*np.pi, np.pi])
    axphi.xaxis.set_ticklabels([r'$-\pi$', r'$-\pi/2.$',0., r'$\pi/2.$',r'$\pi$'])
    f.tight_layout()

    f=plt.figure(); axrho=f.add_subplot(111)
    for ss in shots:
        #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])
        data=dd[ss]['rho']
        axrho.hist(data, histtype='step', color=dd[ss]['lc'], bins=30, lw=2.3, label=dd[ss]['label'])
        #data=dd[ss]['theta']
        #axtheta.hist(data, histtype='step', color=dd[ss]['lc'], bins=20, lw=2.3, label=dd[ss]['label'])

    axrho.legend(loc='upper center')    

    axrho.set_xlabel(r'$\rho$')
    axrho.set_ylabel(r'Lost markers (AU)')
    axrho.grid('on')
    axrho.yaxis.set_ticks([])
    f.tight_layout()

plt.show()
