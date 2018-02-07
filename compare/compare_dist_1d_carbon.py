import ascot_distributions, ascot_particles
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
import scipy.interpolate as interpolate



plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)

# Create your ticker object with M ticks
M = 4
yticks = ticker.MaxNLocator(M)
yticks_m=ticker.MaxNLocator(M*2)
xticks = ticker.MaxNLocator(M)

styles = ['-','--', '-.']
plot_flag = 01
shot = 2
shot = 5
if shot == 2:
    ascot_run = 14
    runs = [7,31]
else:
    ascot_run = 0
    runs = [12,11]

carbon = ["No imp.","Imp."]
power = [0,0,0]
curpow = [0,0,0]


fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(runs[0])+'.h5'
dis = ascot_distributions.distribution_1d(fname)
cur = np.zeros((len(runs), len(dis.volumes)-1), dtype=float)
totcur = np.zeros((len(runs)), dtype=float)
inipow = np.zeros((len(runs)), dtype=float)
endpow = np.zeros((len(runs)), dtype=float)
pow_wall = np.zeros((len(runs)), dtype=float)
totpow_e = np.zeros((len(runs)), dtype=float)
totpow_i = np.zeros((len(runs)), dtype=float)
cur_scal = np.zeros((len(runs)), dtype=float)
powe_scal = np.zeros((len(runs)), dtype=float)
powi_scal = np.zeros((len(runs)), dtype=float)
powc_scal = np.zeros((len(runs)), dtype=float)


pow_i = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
pow_e = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
dens = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
pow_C = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
cur_NNB = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
cur_PNB_par = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
cur_PNB_per = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)

npitch = 10
pitch = np.linspace(-1,1,npitch, dtype=float)
dist_RZE = np.zeros((len(runs), npitch), dtype=float)
tmpdist = np.zeros(20, dtype=float)
for i, el in enumerate(runs):
    fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(el)+'.h5'
    print fname
    dis = ascot_distributions.distribution_1d(fname)
    dis.rhodists()

    #dis.plot_groupcurrent()

    dis.sum_all()
    #cur_NNB[i,:] = dis.data_NNB[:,2]*1e-3
    #cur_PNB_par[i,:] = dis.data_PPAR[:,2]*1e-3
    #cur_PNB_per[i,:] = dis.data_PPERP[:,2]*1e-3
    cur[i,:] = dis.Itot_prof*1e-3
    cur_scal[i] = np.dot(cur[i,:], dis.areas[1:])
    pow_e[i,:] = dis.petot_prof*1e-3
    powe_scal[i] = np.dot(pow_e[i,:], dis.volumes[1:])
    pow_i[i,:] = dis.pitot_prof*1e-3
    powi_scal[i] = np.dot(pow_i[i,:], dis.volumes[1:])
    dens[i,:]  = dis.ntot_prof
    try:
        pow_C[i,:] = dis.pimptot_prof*1e-3
    except:
        pow_C[i,:] = np.zeros(len(dis.volumes)-1, dtype=float)    
    powc_scal[i] = np.dot(pow_C[i,:], dis.volumes[1:])
    

    ptcls = ascot_particles.h5_particles(fname)
    inipow[i] = np.sum(ptcls.data_i['energy']*ptcls.data_i['weight'])*1.6e-19
    endpow[i] = np.sum(ptcls.data_e['energy']*ptcls.data_e['weight'])*1.6e-19
    # find particles that went to the wall
    ind_wall = np.where(ptcls.data_e['endcond'] == 3)
    pow_wall[i] = np.sum(ptcls.data_e['energy'][ind_wall]*ptcls.data_e['weight'][ind_wall])*1.6e-19

#    dis2d = ascot_distributions.frzpe(fname)
#    dis2d.integrate_RZE()
#    if shot ==5 and el == 10:
#        tmpdist = dis2d.f_RZE_int
#    else:
#        dist_RZE[i,:] = dis2d.f_RZE_int
    

for i, el in enumerate(runs):
     print "RUN ", el, " ", carbon[i]
     print "{:05.4f} % absorbed out of {:05.4f} MW".format((inipow[i]-pow_wall[i]-endpow[i])/30e6, 30e6)
     print "{:05.4f} MW lost to the wall after slowing-down".format(pow_wall[i]*1e-6)
     print "{:05.4f} kA NBCD".format(cur_scal[i])
     tmppow = powe_scal[i]+powi_scal[i]+powc_scal[i]
     print "{:05.4f} MW power to e || {:05.4f} % ".format(powe_scal[i]*1e-3,powe_scal[i]/tmppow*100. )
     print "{:05.4f} MW power to i || {:05.4f} % ".format(powi_scal[i]*1e-3,powi_scal[i]/tmppow*100. )
     print "{:05.4f} MW power to C || {:05.4f} % ".format(powc_scal[i]*1e-3,powc_scal[i]/tmppow*100. )

#     print "{:05.4f} kA for injection power".format(totcur[i])
#     print "{:05.4f} MW {:05.4f} % to el. for injection power".format(totpow_e[i],totpow_e[i]/totpow[i]*100.)
#     print "{:05.4f} MW {:05.4f} % to ions for injection power)".format(totpow_i[i], totpow_i[i]/totpow[i]*100.)

if plot_flag == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        ax.plot(dis.rho, -cur[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = str(carbon[i]))
        #ax.plot(dis.rho, cur_NNB[i,:], linewidth = 2.3, linestyle=styles[i],color='r', label = str(carbon[i])+' - I NNB')
        #ax.plot(dis.rho, cur_PNB_par[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = str(carbon[i])+' - I = '+'{:5.2f}'.format(totcur[i])+' kA')

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'j [$kA/m^2$]')
    ax.legend(loc='best')
    #ax.yaxis.set_major_locator(yticks)
    ##ax.yaxis.set_minor_locator(yticks_m)
    #ax.xaxis.set_major_locator(xticks)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        ax.plot(dis.rho, pow_e[i,:]+pow_i[i,:]+pow_C[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = carbon[i]+" Total")

        ax.plot(dis.rho, pow_e[i,:], linewidth = 2.3, linestyle=styles[i],color='r', label = carbon[i]+" e ")
        ax.plot(dis.rho, pow_i[i,:], linewidth = 2.3, linestyle=styles[i],color='g', label = carbon[i]+" D ")
        if shot==5 and el == 11:
            ax.plot(dis.rho, pow_C[i,:], linewidth = 2.3, linestyle=styles[i],color='b', label = carbon[i]+" C ")        

        #ax.plot(dis.rho, pow_i[i,:], linewidth = 2.3, linestyle=styles[i],color='r', label = "ions "+str(en[i])+' keV, run '+str(el))


    ax.yaxis.set_major_locator(yticks)
    #ax.yaxis.set_minor_locator(yticks_m)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'Power density [$kW/m^3$]')
    ax.legend(loc='best')    

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i, el in enumerate(runs):
    #     if el == 10:
    #         ax.plot(np.linspace(-1,1,20), tmpdist, linewidth = 2.3, linestyle=styles[i],color='k', label = str(en[i])+' keV')
    #     else:
    #         ax.plot(pitch, dist_RZE[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = str(en[i])+' keV')
    # ax.yaxis.set_major_locator(yticks)
    # #ax.yaxis.set_minor_locator(yticks_m)
    # ax.xaxis.set_major_locator(xticks)
    # ax.set_xlabel(r'Pitch $v_\parallel/v$')
    # ax.set_ylabel(r'Particle density [$1/m^3$]')
    # ax.legend(loc='best') 
    

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i, el in enumerate(runs):
    #     ax.plot(dis.rho, dens[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = str(carbon[i]))
    # ax.yaxis.set_major_locator(yticks)
    # ax.xaxis.set_major_locator(xticks)
    # ax.set_xlabel(r'$\rho$')
    # ax.set_ylabel(r'Particle density [$1/m^3$]')
    # ax.legend(loc='best')


    plt.show()
