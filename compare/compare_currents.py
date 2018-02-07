import ascot_distributions
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

# Create your ticker object with M ticks
M = 4
yticks = ticker.MaxNLocator(M)
yticks_m=ticker.MaxNLocator(M*2)
xticks = ticker.MaxNLocator(M)

styles = ['-','--', '-.']
plot_flag = 1
shot = 2
ascot_run = 14
runs = [27,28,29]
power = [0,0,0]
curpow = [0,0,0]
en = [500, 350, 250]
for i,e in enumerate(en):
    power[i]=5e6*(e*1000/500000.)**(2.5)
fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(ascot_run)+'.h5'
dis = ascot_distributions.distribution_1d(fname)
cur = np.zeros((len(runs), len(dis.volumes)-1), dtype=float)
totcur = np.zeros((len(runs)), dtype=float)
totpow = np.zeros((len(runs)), dtype=float)
totpow_e = np.zeros((len(runs)), dtype=float)
totpow_i = np.zeros((len(runs)), dtype=float)
pow_i = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)
pow_e = np.zeros((len(runs),len(dis.volumes)-1), dtype=float)

for i, el in enumerate(runs):
    fname = '/home/vallar/ASCOT/runs/JT60SA/'+"{:03d}".format(shot)+'/ascot_'+"{:03d}".format(shot)+"{:03d}".format(el)+'.h5'
    dis = ascot_distributions.distribution_1d(fname)
    dis.rhodists()
    dis.group_beams()
    cur[i,:] = dis.data_NNB[:,2]*1e-3
    pow_e[i,:] = dis.data_NNB[:,21]*1e-6
    pow_i[i,:] = dis.data_NNB[:,22]*1e-6
    dis.calculate_scalar()
    totcur[i] = dis.I_tot*1e-3
    totpow_e[i] = (dis.pe_beams[-2]+dis.pe_beams[-1])*1e-6
    totpow_i[i] = (dis.pi_beams[-2]+dis.pi_beams[-1])*1e-6
    totpow[i] = totpow_e[i]+totpow_i[i]
    curpow[i] = totcur[i]/(power[i]*1e-6)
    
for i, el in enumerate(runs):
    print "RUN ", el
    print "{:05.4f} kA for injection power {:05.4f} ({:05.4f} in total)".format(totcur[i], power[i]*1e-6, 2*power[i]*1e-6)
    print "{:05.4f} kA/MW for injection power {:05.4f} ({:05.4f} in total)".format(curpow[i],  power[i]*1e-6, 2*power[i]*1e-6)
    print "{:05.4f} MW {:05.4f} % to el. for injection power {:05.4f} ({:05.4f} in total)".format(totpow_e[i],totpow_e[i]/totpow[i]*100. , power[i]*1e-6, 2*power[i]*1e-6)
    print "{:05.4f} MW {:05.4f} % to ions for injection power {:05.4f} ({:05.4f} in total)".format(totpow_i[i], totpow_i[i]/totpow[i]*100., power[i]*1e-6, 2*power[i]*1e-6)
    print "{:05.4f} MW for injection power {:05.4f} ({:05.4f} in total)".format(totpow[i],power[i]*1e-6, 2*power[i]*1e-6)

if plot_flag == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        ax.plot(dis.rho, cur[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = str(en[i])+' keV - I = '+'{:5.2f}'.format(totcur[i])+' kA')
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'j [$kA/m^2$]')
    ax.legend(loc='best')
    #ax.yaxis.set_major_locator(yticks)
    ##ax.yaxis.set_minor_locator(yticks_m)
    #ax.xaxis.set_major_locator(xticks)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, el in enumerate(runs):
        ax.plot(dis.rho, pow_e[i,:], linewidth = 2.3, linestyle=styles[i],color='r', label = "Electrons "+str(en[i])+' keV')
        #ax.plot(dis.rho, pow_i[i,:], linewidth = 2.3, linestyle=styles[i],color='r', label = "ions "+str(en[i])+' keV, run '+str(el))
        ax.plot(dis.rho, pow_e[i,:]+pow_i[i,:], linewidth = 2.3, linestyle=styles[i],color='k', label = "Total "+str(en[i])+' keV')


    ax.yaxis.set_major_locator(yticks)
    #ax.yaxis.set_minor_locator(yticks_m)
    ax.xaxis.set_major_locator(xticks)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'Power density [$MW/m^3$]')
    ax.legend(loc='best')    
    plt.show()
