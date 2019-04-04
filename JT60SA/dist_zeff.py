"""
Script to join different profiles (w 1e6 particles for each set of beams)
"""
import a4py.classes.distributions as ascot_distributions
import a4py.classes.particles as ascot_particles
import a4py.classes.distributions_2d as d2d
import a4py.classes.orbits as orbits

from utils.plot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

def interp1d(x,y,xnew):
    p = interp.interp1d(x,y)
    ynew = p(xnew)
    return ynew

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
    dd[el]['dE'] = p.data_e['energy'][ind]
    dd[el]['particles'] = p
    p.endcondition()
    #dist2d = d2d.frzpe(folder+fname)
    #dd[el]['d2d'] = dist2d
    dist = ascot_distributions.SA_1d(folder+fname)
    dist._group_beams()
    dd[el]['rho'] = dist.rho
    dd[el]['j'] = dist.slices_summed[0,:,12]
    dd[el]['pel'] = dist.slices_summed[0,:,dist.peind]
    dd[el]['pi'] = dist.slices_summed[0,:,dist.peind+1]
    dd[el]['pr'] = dist.slices_summed[0,:,13]
    dd[el]['n'] = dist.slices_summed[0,:,0]
    dd[el]['pi2'] = dist.slices_summed[0,:,dist.peind+2]
    dist.calc_eff();
    dd[el]['eff'] = dist.eff
    #_,_,_,dd[el]['taus'] = dist._ecrit(E0=dd[el]['E'])
    
    #dist.print_scalars()
    dd[el]['gi'] = 1-dist.pe/(dist.pe+dist.pi1+dist.pi2)
    dist.print_scalars()
    _,ec,_, taus = dist._ecrit()
    print()
    print('slowing-down time for 500 keV: ', taus, ' Ec= ', ec)
    _,ec,_, taus = dist._ecrit(E0=85000)
    print('slowing-down time for 85 keV: ', taus, ' Ec= ', ec)
    print()
f=plt.figure(); axetime=f.add_subplot(111)
bins = np.linspace(10, 90, 20)
for ss in shots:
    #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])

    data=dd[ss]['dE']*1e-3
    axetime.hist(data, histtype='step', color=dd[ss]['lc'], bins=20, lw=2.3, label=dd[ss]['label'])
axetime.legend(loc='upper center')
axetime.set_xlabel(r'E$_{final}$ [kev]')
axetime.set_ylabel(r'Makers at the wall')
axetime.grid('on')
f.tight_layout()
if False:
    f=plt.figure(); axpitch=f.add_subplot(111)
    bins = np.linspace(10, 90, 20)
    for ss in shots:
        #_plot_1d(dd[ss]['dE'], hist=1, ax=axetime, color=dd[ss]['lc'])
        dd[ss]['d2d'].plot_spaceE(ax=axpitch, color=dd[ss]['lc'], label=ss)
    axpitch.legend(loc='lower center')
    axpitch.set_xlabel(r'$\xi$')
    axpitch.set_ylabel(r'Normalized fast ion density', fontsize='xx-large')
    axpitch.yaxis.set_ticks([])
    axpitch.grid('on')
    f.tight_layout()


if True:
    rho = dd['Ref.']['rho']
    data_labels = ['Ref.', 'Acc.','Edge']
    ylab = [r'j (kA/$m^2$)', r'$P_{el}\,(kW/m^3)$', r'$P_{D}\,(kW/m^3)$', r'p (kN/$m^2$)', r'n (1/$m^3$)', r'$P_{C}\,(kW/m^3)$']
    factor = [-1e-3, 1e-3, 1e-3, 1e-3, 1, 1e-3]
    col= [dd['Ref.']['lc'],dd['Core']['lc'],dd['Edge']['lc']]
    for i, el in enumerate(['j', 'pel', 'pi', 'pr','n', 'pi2']):
        #tot = factor[i]*dd['500'][el].T+factor[i]*dd['350'][el].T+factor[i]*dd['250'][el].T
        data = np.array([rho, factor[i]*dd['Ref.'][el].T, factor[i]*dd['Core'][el].T, factor[i]*dd['Edge'][el].T])
        if el=='pel' or el=='pi' or el == 'pi2':
            if shot==4 or shot==5:
                ylim=[0,350]
            else:
                ylim=[0,200]
        else:
            ylim=0
        plot_article(3,data, data_labels=data_labels,xlabel=r'$\mathbf{\rho}$', ylabel=ylab[i], ylim=ylim, col=col)

plt.show()
