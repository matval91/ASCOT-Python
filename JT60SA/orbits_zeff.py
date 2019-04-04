import a4py.classes.orbits as orbits

from utils.plot_utils import _plot_1d, plot_article, limit_labels
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import pickle

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
    dd['Ref.']['fname'] = '/orbits_ref.dat'
    dd['Core']['fname'] = '_core/orbits_core.dat'
    dd['Edge']['fname'] = '_edge/orbits_edge.dat'
    dd['Ref.']['fname_surf'] = '/ascot_dist.h5'
    dd['Core']['fname_surf'] = '_core/ascot_dist.h5'
    dd['Edge']['fname_surf'] = '_edge/ascot_dist.h5'

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
for el in dd:
    fname = dd[el]['fname']
    fname_surf = dd[el]['fname_surf']

    o = orbits.SA_orbits(folder+fname, fname_surf=folder+fname_surf) 
    dd[el]['o'] = o

figr=plt.figure(); axrho=figr.add_subplot(111)
figp=plt.figure(); axpitch=figp.add_subplot(111)
figw=plt.figure(); axw=figw.add_subplot(211); axw_n = figw.add_subplot(212)
lim_banana_nnb=7.5*0.01
for el in dd:
    o = dd[el]['o']
    rho = o.data_i['rho'][np.where(o.trappind==1)[0]]
    axrho.hist(rho, histtype='step', color=dd[el]['lc'], label=el, lw=2.3, density=True, bins=40)
    pitch = o.data_i['pitch'][np.where(o.trappind==1)[0]]
    axpitch.hist(pitch, histtype='step', color=dd[el]['lc'],  label=el, lw=2.3, density=True, bins=40)
    indp = np.where(o.w_banana<lim_banana_nnb)[0]
    axw.hist(o.w_banana[indp]*100., histtype='step', label=el, color=dd[el]['lc'], lw=2.3, density=True, bins=40)
    indn = np.where(o.w_banana>lim_banana_nnb)[0]
    axw_n.hist(o.w_banana[indn]*100., histtype='step', label=el, color=dd[el]['lc'], lw=2.3, density=True, bins=40)

limit_labels(axw, r'$w_{banana}$ [cm]'); limit_labels(axw_n, r'$w_{banana}$ [cm]')  
limit_labels(axrho, r'$\rho$')  
limit_labels(axpitch, r'$\xi$')      
figr.tight_layout()
figp.tight_layout()
plt.show()

dddd={key:[] for key in dd}
for el in dd:
    o = dd[el]['o']
    rho = o.data_i['rho'][np.where(o.trappind==1)[0]]
    pitch = o.data_i['pitch'][np.where(o.trappind==1)[0]]
    wbanana = o.w_banana
    
    dddd[el].append({'rho':rho, 'pitch':pitch, 'wbanana':wbanana})
with open('data.txt', 'w') as outfile: 
    json.dumps(dddd, outfile, cls=NumpyEncoder)

output = open('myfile.pkl', 'wb')
pickle.dump(dddd, output)
output.close()
