import a5py.ascot5io.ascot5 as a5
from utils.plot_utils import _plot_2d

f=a5.Ascot('ascot.h5')
e=aa.active.endstate.read()
endcond=e['endcond']
ind=np.where(e['endcond']==8)[0] #wall collision endstate
R=e['R'][ind]; z=e['z'][ind]
v_end = np.sqrt(e['vR']**2+e['vz']**2+e['vphi']**2)
e_end = e['mass']*1.66054e-27*0.5*v_end*v_end/1.602e-19

i=aa.active.inistate.read()
rho = i['rho']
v_ini = np.sqrt(i['vR']**2+i['vz']**2+i['vphi']**2)
pitch = i['vpar']/v_ini
_plot_2d(rho[ind], pitch[ind], xlabel=r'$\rho$', ylabel=r'$\xi$', scatter=1)

_plot_2d(e_end[ind]*1e-3,e['time'][ind], xlabel=r'$E_{end} [keV]$', ylabel=r't [s]', scatter=1)
