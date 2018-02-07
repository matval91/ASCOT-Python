import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, colors
import collections
import math
import h5py
from collections import OrderedDict
import ascot_distributions

cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)

class particles:
    """
    SUPERCLASS
    
    """
    def __init__(self):
        
        # Initialisation
        self.npart   = 0
        self.nfields = 0
        self.field   = []
        self.unit    = []
        self.data_i  = collections.defaultdict(list)
        self.data_e  = collections.defaultdict(list)
        self.R_w = []
        self.z_w = []

    def _calc_binarea(self,x,y):
        """
        hidden method to calculate the XY area of the bin for XY ionisation plot
        """
        Dx = x[1]-x[0]
        Dy = y[1]-y[0]
        area = Dx*Dy
        return area

    def _calc_weight(self):
        """
        hidden method to calculate the weight of the particles for each origin: in the inistate,
        in the final state and the total one.
        """
        try:
            self.data_i['weight'][:].mean()
        except:
            print "No weights available"
            return
        
        self._calc_originbeam()
        self.w_i = np.zeros((len(self.origins)), dtype=float)
        self.w_e = np.zeros((len(self.origins)), dtype=float)
        self.w   = np.zeros((len(self.origins)), dtype=float)

        for ind_i, i in enumerate(self.origins):
            ind = self.data_i['origin'][:]==i
            self.w_i[ind_i] = sum(self.data_i['weight'][ind])
            ind = self.data_e['origin'][:]==i
            self.w_e[ind_i] = sum(self.data_e['weight'][ind])  
            self.w[ind_i]  = self.w_i[ind_i]+self.w_e[ind_i]
            
    def plot_RZpart(self):
        """
        Method to plot R vs z of the ionised particles, useful mostly with bbnbi
        """
        x=self.data_i['Rprt']
        y=self.data_i['zprt']

        f=plt.figure()
        ax = f.add_subplot(111)
        #hb = ax.hexbin(x, y, gridsize=100, cmap=my_cmap)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap)
        cb = f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('R')
        ax.set_ylabel('z')

        if len(self.R_w)==0:
            self._readwall()
        ax.plot(self.R_w , self.z_w, 'k', linewidth=3)

        ax.axis('equal')
        ax.axis([min(self.R_w), max(self.R_w), min(self.z_w), max(self.z_w)])
        #ax.set_xrange([min(self.R_w), max(self.R_w)])
        #ax.set_yrange([min(self.z_w), max(self.z_w)])                   
        plt.show()


    def plot_XYpart_singlebeam(self):
        """
        Method to plot XY of ionisation, without difference between the beams
        """

        x=np.zeros(self.npart)
        y=np.zeros(self.npart)
        
        for i, el in enumerate(self.data_i['Rprt']):
            x[i]=el*math.cos(self.data_i['phiprt'][i])
            y[i]=el*math.sin(self.data_i['phiprt'][i])
            
        f=plt.figure()
        ax=f.add_subplot(111)
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap, normed=True)
        cb = f.colorbar(hb[3], ax=ax)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        theta=np.arange(0,6.3,0.02*6.28)
        if len(self.R_w)==0:
            self._readwall()
        ax.plot(np.min(self.R_w)*np.cos(theta) , np.min(self.R_w)*np.sin(theta), 'k', linewidth=3)
        ax.plot(np.max(self.R_w)*np.cos(theta) , np.max(self.R_w)*np.sin(theta), 'k', linewidth=3)
        plt.show()

    def plot_ioniz_beams(self):
        """
        Method to plot XY of ionisation FOR EACH BEAM
        """
        try:
            self.beamorigindict.keys()
        except:
            self._calc_originbeam()

            
        #x_range=[-5,5]
        axxy = plt.subplot2grid((1,2),(0,0))
        axrz = plt.subplot2grid((1,2),(0,1))
        ind = np.zeros((2,self.npart), dtype=bool)
        
        for i,el in enumerate(self.beamorigindict.keys()):
            if el != 'NNBI_U' and el!='NNBI_L':
                ind[0,:] = self.data_i['origin']==self.beamorigindict[el][0]
                ind[1,:] = self.data_i['origin']==self.beamorigindict[el][1]
            elif el=='NNBI_L':
                continue
            else:
                ind[0,:] = self.data_i['origin']==self.beamorigindict[el]
                ind[1,:] = self.data_i['origin']==self.beamorigindict['NNBI_L']


            for jj in (0,1):
                R   = self.data_i['Rprt'][ind[jj,:]]
                ang = self.data_i['phiprt'][ind[jj,:]]
                z   = self.data_i['zprt'][ind[jj,:]]
                x = np.zeros(len(R))
                y = np.zeros(len(R))
                for j, el in enumerate(R):
                    x[j] = el*math.cos(ang[j])
                    y[j] = el*math.sin(ang[j])
                axxy.scatter(x,y,c=col[i])            
                axrz.scatter(R,z,c=col[i])
        
##         axxy.('X')
##         axxy.ylabel('Y')
##         axrz.xlabel('R')
##         axrz.ylabel('z')
        theta=np.arange(0,6.29,0.02*6.28)
        if len(self.R_w)==0:
            self._readwall()
        axxy.plot(np.min(self.R_w)*np.cos(theta) , np.min(self.R_w)*np.sin(theta), 'm')
        axxy.plot(np.max(self.R_w)*np.cos(theta) , np.max(self.R_w)*np.sin(theta), 'm')
        axrz.plot(self.R_w , self.z_w, 'm')
        axxy.axis([-5,5,-5,5])
        axrz.axis([0,5,-4,4])
        axxy.axis('equal')
        axrz.axis('equal')
        
        plt.show()


    def losses_prediction(self):
        endstate=self.data_e['endcond']
        R   = self.data_i['R']
        phi = self.data_i['phi']
        plt.figure()
        plt.scatter(R*np.cos(phi), R*np.sin(phi))
        ind_th   = endstate == 4
        ind_wall = endstate == 3
        ind_emin = endstate == 2
        plt.figure()
        plt.scatter(self.data_i['pitch'][ind_wall], self.data_i['energy'][ind_wall],c='blue', label='wall', marker='x')
        plt.scatter(self.data_i['pitch'][ind_th],   self.data_i['energy'][ind_th], c='red',  label='th', marker='x')
        plt.scatter(self.data_i['pitch'][ind_emin], self.data_i['energy'][ind_emin],c='green', label='emin', marker='x')
        plt.legend(loc='best')

        plt.figure()
        plt.scatter(R[ind_wall]*np.cos(phi[ind_wall]), R[ind_wall]*np.sin(phi[ind_wall]),c='blue', label='wall', marker='x')
        plt.scatter(R[ind_th]*np.cos(phi[ind_th]), R[ind_th]*np.sin(phi[ind_th]),c='red', label='th', marker='x')
        plt.scatter(R[ind_emin]*np.cos(phi[ind_emin]), R[ind_emin]*np.sin(phi[ind_emin]),c='green', label='emin', marker='x')
        
        plt.legend(loc='best')


    def _readwall(self):
        """
        Hidden method to read the wall
        """
        in_w_fname = 'input.wall_2d'
        try:
            wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=1)
        except:
            in_w_fname = '../input.wall_2d'
            wall = np.loadtxt( in_w_fname, dtype=float, unpack=True, skiprows=1)
        self.R_w = wall[0,:]
        self.z_w = wall[1,:]
        self.R_w = np.array(self.R_w)
        self.z_w = np.array(self.z_w)

    def _calc_originbeam(self):
        """
        Method to translate from origin in bbnbi.h5 file to beam identification
        Now only with JT60-SA (and still to recognize A and B units)
        """
        self.beamorigindict=OrderedDict()
        self.beamorigindict={'1':[45,46]    , '2':[47,48]     ,'3':[133,134]   ,'4':[135,136],\
                             '5':[221,222]  , '6':[223,224]   ,'7':[309,310]   ,'8':[311,312],\
                             '9':[3637,3638], '10':[3639,3640],'13':[5253,5254],'14':[5255,5256],\
                             'NNBI_U':[3031], 'NNBI_L':[3032]}
        
        self.origindict={'45':'1_1', '46':'1_2', '47':'2_1', '48':'2_2',\
                         '133':'3_1', '134':'3_2', '135':'4_1', '136':'4_2',\
                         '221':'5_1', '222':'5_2', '223':'6_1', '224':'6_2',\
                         '309':'7_1', '310':'7_2', '311':'8_1', '312':'8_2',\
                         '3637':'9_1', '3638':'9_2', '3639':'10_1', '3640':'10_2',\
                         '5253':'13_1', '5254':'13_2', '5255':'14_1', '5256':'14_2',\
                         '3031':'NNB_1', '3032':'NNB_2'}
        self.beamlabel=['1','2','3','4','5','6','7','8','9','10',\
                        '13','14','NNB_U','NNB_L']
        self.beamlabel_num=['1','2','3','4','5','6','7','8','9','10',\
                            '13','14','99','101']
        self.beamorigin = np.array([45,46,47,48,133,134,135,136,\
                           221,222,223,224,309,310,311,312,\
                           3637,3638,3639,3640,5253,5254,5255,5256,\
                           3031,3032], dtype=int)
    

        self.originpower={'45':1, '46':1, '47':1, '48':1,\
                         '133':1, '134':1, '135':1, '136':1,\
                         '221':1, '222':1, '223':1, '224':1,\
                         '309':1, '310':1, '311':1, '312':1,\
                         '3637':1, '3638':1, '3639':1, '3640':1,\
                         '5253':1, '5254':1, '5255':1, '5256':1,\
                         '3031':5, '3032':5}


        #orderedorigins is the array where you can retrieve the correct beam used in the distributions
        # so in orderedorigins we will find e.g. [13,7,9,11,15,4,8,9,...]
        # and you get data of beam 1 using orderedorigins[0], etc.
##         self.orderedorigins=np.zeros(len(self._h5origins), dtype=int)
##         for i,el in enumerate(self.beamorigin):
##             ind = np.where(self._h5origins==el)
##             self.orderedorigins[i] = int(ind[0])


    def endcondition(self):
        """
        Computes the endcondition of the particles and stores a dict with the amount of MC particles with that final condition
        """
        errendcond  = {'aborted':-2, 'rejected':-1}
        counter = 0
        for key in errendcond:
            ind = np.where(self.data_e['endcond']==errendcond[key])[0]
            if len(ind)!=0:
                counter += len(ind)
                print("STRANGE END CONDITION! {} ({:d}) : {:d} particles".format(key, errcond[key], len(ind)))

        physcond = {'none':0,'tmax':1,'emin':2,'wall':3,'thermalisation':4}
        for key in physcond:
            ind = np.where(self.data_e['endcond']== physcond[key])[0]
            counter += len(ind)
            print("{} ({:2d}) : {:4d} particles".format(key, physcond[key], len(ind)))

        infcond = {'cputmax':17, 'outsidemaxrho':10, 'outsidewall':16}
        for key in infcond:
            ind = np.where(self.data_e['endcond']== infcond[key])[0]
            counter += len(ind)
            print("{} ({:2d}) : {:4d} particles".format(key, infcond[key], len(ind)))    

        print "Total particles counted :", counter/len(self.data_e['endcond'])*100., ' %'
        

class dat_particles(particles):
    """
    Class (inherited from particles) handling dat files (e.g. input.particles, test_ascot.orbits.dat)
    """
    def __init__(self, fname):
        """
        READ ASCII FILE
        """
        particles.__init__(self)
        self.file = fname
        in_f  = open(self.file)
        lines = in_f.readlines()[3:]
        in_f.close()
        
        #Read lines of comment
        n_comm = int(lines[0].split()[0])
        
        # read number of particles
        self.npart = int(lines[n_comm+2].split()[0])

        # Number of fields
        self.nfields = int(lines[n_comm+4].split()[0])
        line_fields  = n_comm+5
        # Read and store the fields and their unit
        for i in range(self.nfields):
            self.field.append(lines[line_fields+i].split()[0])
            self.unit.append(lines[line_fields+i].split()[-1][1:-1])

        # Store actual data 
        part = lines[line_fields+self.nfields+1:-1]
        for rowt in part:
            row = rowt.split()
            for i, el in enumerate(row):
                self.data_i[self.field[i]].append(el)

        for key in self.field:
            if key=='id':
                self.data_i[key] = np.array(self.data_i[key], dtype=int)
            else:
                self.data_i[key] = np.array(self.data_i[key], dtype=float)

        #CONVERT PHI FROM DEG TO RAD
        self.data_i['phiprt']=self.data_i['phiprt']*math.pi/180.
        self.origins = np.array(list(set(self.data_i['origin'])))
        print "FIELDS IN FILE ",self.file," :"
        print self.field
        
        if self.npart==-1:
            self.npart = np.max(self.data_i['id'])


class h5_particles(particles):
    """
    Class (inherited from particles) handling h5 files (e.g. bbnbi.h5, ascot.h5)
    """
    def __init__(self, fname):
        """
        READ *.h5 FILE
        """
        particles.__init__(self)
        self.file = fname
        dataf = h5py.File(self.file)
        try:
            self.endgroup = 'shinethr' #this is for bbnbi
            self.field = dataf[self.endgroup].keys()
        except:
            self.endgroup = 'endstate' #this is for ascot
            self.field = dataf[self.endgroup].keys()
            
        self.nfields=len(self.field)
        for key in self.field:
            if key=='id':
                self.data_i[key] = np.array(dataf['inistate/'+key], dtype=int)
                self.data_e[key] = np.array(dataf[self.endgroup+'/'+key], dtype=int)
            else:
                self.data_i[key] = np.array(dataf["inistate/"+key], dtype=float)
                self.data_e[key] = np.array(dataf[self.endgroup+"/"+key], dtype=float)

        if self.endgroup == 'endstate':
            self._h5origins = dataf['species/testParticle/origin'].value
        elif self.endgroup == 'shinethr':
            self.origins = np.sort(np.array(list(set(self.data_i['origin']))))

        # evaluate number of particles
        self.npart = np.max(self.data_i['id'])

        #CONVERT PHI FROM DEG TO RAD
        self.data_i['phiprt']=self.data_i['phiprt']*math.pi/180.
        self.data_e['phiprt']=self.data_e['phiprt']*math.pi/180.

        print "FIELDS IN FILE ",self.file,", section inistate and ",self.endgroup," :"
        print self.field


    def TCV_calc_shinethr(self):
        """
        Method to calculate the shine-through with the weights
        """
        
        if self.endgroup != 'shinethr':
            print "WARNING! Check input file, which is ", self.file

        self._calc_weight()
        self.shinethr=np.zeros((1), dtype=float)
        self.shinethr_abs=np.zeros((1), dtype=float)
        power = 1e6

        e = self.data_e['energy']
        e = e*1.612e-19
        w = self.data_e['weight']
        self.shinethr_abs= np.sum(e*w)
        self.shinethr=np.sum(e*w)/power



        print "TCV Shine-through:", "{0:5f} %".format(self.shinethr*100),\
                  "||  {0:3f} W".format(self.shinethr_abs)



    def calc_shinethr(self):
        """
        Method to calculate the shine-through with the weights
        """
        
        if self.endgroup != 'shinethr':
            print "WARNING! Check input file, which is ", self.file

        self._calc_weight()
        self.shinethr=np.zeros((len(self.origins)), dtype=float)
        self.shinethr_abs=np.zeros((len(self.origins)), dtype=float)

        for ind_i,i in enumerate(self.origins):
            ind = self.data_e['origin'][:]==i
            if len(self.data_e['origin'][ind])==0:
                w=0
                e=0
                power = 1
            else:   
                e = self.data_e['energy'][ind][0]
                e = e*1.612e-19
                w = np.sum(self.data_e['weight'][ind_i])
                #wtot = self.w[ind_i]
                power = self.originpower[str(int(i))]*1e6
            self.shinethr[ind_i]=float(e)*w/power
            self.shinethr_abs[ind_i] = float(e)*w
            print "NBI:","{}".format(self.origindict[str(int(i))]),\
                  "Shine-through:", "{0:5f} %".format(self.shinethr[ind_i]*100),\
                  "||  {0:3f} W".format(e*w)


    def particle_status(self):
        """
        With this function you can plot the histogram with thermalised and lost particles
        """
        # text font for plot
        plt.rc('font', family='serif', serif='Palatino')
        #plt.rc('text', usetex=True)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=20)
        #=========================

        
        endstate=self.data_e['endcond']
        self._calc_originbeam()
        self.npart_arr = np.zeros(14, dtype=float)
        self.endstatus=np.zeros((14, 5), dtype=float)
        order=[]
        frac_wall = np.zeros(14, dtype=float)
        frac_rho = np.zeros(14, dtype=float)
        frac_tot_rho, frac_tot_wall = 0,0
        keys=['1','2','3','4','5','6','7','8','9','10','13','14','NNBI_U','NNBI_L']
        for ind, el in enumerate(keys):
            n_part=0
            n_wall=0
            n_th=0
            n_emin=0
            n_rho = 0
            order.append('$'+str(el)+'$')
            for jj in range(self.npart):

                if self.data_e['origin'][jj] in self.beamorigindict[el]:
                    n_part+=1
                    if endstate[jj] == 3: #WALL
                        n_wall+=1
                    elif endstate[jj] == 4: #THERMALISATION
                        n_th+=1
                    elif endstate[jj] == 2: # EMIN
                        n_emin+=1
                    elif endstate[jj] == 10: # OUTSIDE MAX RHO
                        n_rho+=1
            self.npart_arr[ind] = n_part
            self.endstatus[ind, 0] = n_part-n_wall-n_th-n_emin-n_rho
            self.endstatus[ind, 1] = n_wall
            self.endstatus[ind, 2] = n_th
            self.endstatus[ind, 3] = n_emin   
            self.endstatus[ind, 4] = n_rho                                                                                        



        xlab=['1','2','3','4','5','6','7','8','9','10','13','14','N-U','N-L']
        width=0.8

        fig=plt.figure(figsize=(12,5))
        ax=fig.add_subplot(111)
        ax.bar(np.arange(14)+0.5,self.endstatus[:,2], width, color='g', label='Th.')
        if n_wall!=0:      
            ax.bar(np.arange(14)+0.5,self.endstatus[:,1], width, bottom=self.endstatus[:,2], color='r', label='Wall')
            frac_wall = self.endstatus[:,1]/self.npart_arr[:]
            frac_tot_wall = np.sum(self.endstatus[:,1])/np.sum(self.npart_arr)
        ax.bar(np.arange(14)+0.5,self.endstatus[:,3], width, bottom=self.endstatus[:,2]+self.endstatus[:,1], color='y', label='$E_{min}$')        
        if n_rho!=0:        
            ax.bar(np.arange(14)+0.5,self.endstatus[:,4], width, bottom=self.endstatus[:,3]+self.endstatus[:,2]+self.endstatus[:,1], color='b', label=r'$\rho_{MAX}$')        
            frac_rho = self.endstatus[:,4]/self.npart_arr[:]
            frac_tot_rho = np.sum(self.endstatus[:,4])/np.sum(self.npart_arr)

        #ax.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14], xlab)
        # Create your ticker object with M ticks
        M = 4
        yticks = ticker.MaxNLocator(M)
        # Set the yaxis major locator using your ticker object. You can also choose the minor
        # tick positions with set_minor_locator.
        ax.yaxis.set_major_locator(yticks)

        #SET TICKS LABELS ON X AXIS
        ax.xaxis.set_ticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
        ax.xaxis.set_ticklabels(xlab)

        
        plt.ylabel('Particles')
        plt.xlabel('Beam index')
        plt.legend(loc='best', prop={'size':20})

        print "Fraction to the wall:",frac_wall
        print "Fraction outside the rho:", frac_rho
        print "Total fraction to wall:", frac_tot_wall
        print "Total fraction outside rho:", frac_tot_rho
