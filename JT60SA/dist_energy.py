"""
Script to join different profiles (w 1e6 particles for each set of beams)
"""
import ascot_distributions, ascot_particles
from ascot_utils import _plot_1d, plot_article
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.interpolate as interp


plt.rc('font', weight='bold')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
plt.rc('figure', facecolor='white')
plt.rc('legend', fontsize=20)


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
elif shot==4:
    folder += '004/'
    dd = {
        '220': {
            'fname':'ascot_004017.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':220e3},
        '400': {
            'fname':'ascot_004016.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':400e3},
        '500': {
            'fname':'ascot_004015.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':500e3}
    }
elif shot==5:
    folder += '005/'
    dd = {
        '500': {
            'fname':'ascot_005067.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':500e3},
        '350': {
            'fname':'ascot_005075.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':350e3},
        '250': {
            'fname':'ascot_005073.h5',
            'j': [],
            'pel':[],
            'pi':[],     
            'pr':[], 'n':[],           
            'rho':[], 'gi':0,
            'eff':0, 'taus':0,
            'E':250e3}
    }

#f=plt.figure()
#a=f.add_subplot(111)
for el in dd:
    fname = dd[el]['fname']
    dist = ascot_distributions.SA_1d(folder+fname)
    p=ascot_particles.SA_iniend(folder+fname)
    p.endcondition()
    dist.group_beams()
    dd[el]['rho'] = dist.rho
    dd[el]['j'] = dist.data_NNB[:,12]
    dd[el]['pel'] = dist.data_NNB[:,dist.peind-1]
    dd[el]['pi'] = dist.data_NNB[:,dist.peind]
    dd[el]['pr'] = dist.data_NNB[:,13]
    dd[el]['n'] = dist.data_NNB[:,0]

    dist.calc_eff();
    dd[el]['eff'] = dist.eff*10.
    _,_,_,dd[el]['taus'] = dist._ecrit(E0=dd[el]['E'])
    
    dist.print_scalars()
    dd[el]['gi'] = (1-dist.pe/(dist.pe+dist.pi1+dist.pi2))*100.

    #d2d=ascot_distributions.frzpe(folder+fname)
    #d2d.plot_spaceE(ax=a, label=el)


#a.legend(loc='best')

rho = dd['500']['rho']


data_labels = ['500','350','250']
ylab = [r'j (kA/$m^2$)', r'$P_{el}\,(kW/m^3)$', r'$P_{i}\,(kW/m^3)$', r'p (kN/$m^2$)', r'n (1/$m^3$)']
ylab = [r'j (kA/$m^2$)', r'p (kN/$m^2$)']
factor = [-1e-3, 1e-3, 1e-3, 1e-3, 1]
for i, el in enumerate(['j','pr']):
    if shot==4:
        tot = factor[i]*dd['500'][el].T+factor[i]*dd['400'][el].T+factor[i]*dd['220'][el].T
        data = np.array([rho, tot, factor[i]*dd['500'][el].T, factor[i]*dd['400'][el].T, factor[i]*dd['220'][el].T])
    else:
        tot = factor[i]*dd['500'][el].T+factor[i]*dd['350'][el].T+factor[i]*dd['250'][el].T
        data = np.array([rho,factor[i]*dd['500'][el].T, factor[i]*dd['350'][el].T, factor[i]*dd['250'][el].T])
    if el=='pel' or el=='pi':
        if shot==4 or shot==5:
            ylim=[0,300]
        else:
            ylim=[0,200]
    else:
        ylim=0
    plot_article(3,data, data_labels=data_labels,xlabel=r'$\mathbf{\rho}$', ylabel=ylab[i], ylim=ylim)

if 1>2:
	f=plt.figure(); ax=f.add_subplot(211); ax2=f.add_subplot(212)
        if shot==4:
            E=[220,400,500]
            eff = [dd['220']['eff'],dd['400']['eff'],dd['500']['eff']]
            gi = [dd['220']['gi'],dd['400']['gi'],dd['500']['gi']]
            taus = [dd['220']['taus'],dd['400']['taus'],dd['500']['taus']]
        else:
            E=[250,350,500]
            eff = [dd['250']['eff'],dd['350']['eff'],dd['500']['eff']]
            gi = [dd['250']['gi'],dd['350']['gi'],dd['500']['gi']]
            taus = [dd['250']['taus'],dd['350']['taus'],dd['500']['taus']]

	ax.scatter(E, gi, 50, color='k')
	ax.set_ylim([30,60]); ax.set_yticks([30,40,50,60])
	ax.set_xticks([250,300,350,400,450,500])
	ax.set_xlim([220, 520])
	ax.set_ylabel(r'$G_i (\%)$')
	ax.grid('on')
	ax.set_xticklabels([])
	#Removing first point of y-axis
	#plt.setp(ax.get_yticklabels()[0], visible=False) 

	axb = ax.twinx()
	axb.scatter(E, eff, 50, color='r')
	axb.set_ylabel(r'$\eta [10^{19} A/Wm^2]$', color='r')
	axb.tick_params('y', colors='r')
	axb.set_yticks([0.5,0.7,0.9,1.1])
	#plt.setp(axb.get_yticklabels()[0], visible=False) 
	axb.set_xlim([200, 520])

	ax2.scatter(E, taus,50, color='k')
	ax2.set_xticks([250,300,350,400,450,500])
	ax2.set_xlim([200, 520])
	ax2.set_yticks([0.2, 0.25, 0.3, 0.35])
        ax2.set_ylim([0.2, 0.35])
	ax2.set_xlabel(r'E [keV]')
	ax2.set_ylabel(r'$\tau_s$ [s]')
	ax2.grid('on')
	#Removing first point of y-axis
	#plt.setp(ax2.get_yticklabels()[0], visible=False) 
	f.tight_layout()
plt.show()
