import matplotlib.pyplot as plt
from matplotlib import ticker, colors
import numpy as np
colours_old = ['k', 'g', 'b', 'r', 'c']
colours = ['k', 'r', 'b', 'g', 'c']

styles = ['-','--','-.']

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

dpi=800

def common_style():
    """
    Defining common styles using plt.rc
    """
    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
    plt.rc('figure', facecolor='white')
    plt.rc('legend', fontsize=20)

def limit_labels(ax, xlabel='', ylabel='', title='', M=5):
    """
    Limiting labels of the axis to 4 elements and setting grid
    """
    #==============================================
    # SET TICK LOCATION
    #==============================================
    
    # Create your ticker object with M ticks
    yticks = ticker.MaxNLocator(M)
    xticks = ticker.MaxNLocator(M)
    # tick positions with set_minor_locator.
    ax.yaxis.set_major_locator(yticks)
    #ax.yaxis.set_minor_locator(yticks_m)
    ax.xaxis.set_major_locator(xticks)
    #==============================================
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid('on', alpha=0.6)
    #Removing first point of y-axis
    plt.setp(ax.get_yticklabels()[0], visible=False) 


def plot_article(n_lines, data, data_labels, xlabel, ylabel, title='', ax=0, ylim=0, fname='', col=['']):
    """
    """
    common_style()
    #===============================
    figsize=[10,8]; flag_label=0
    n_l_oplot=0
    if ax==0:
        fig=plt.figure(figsize=figsize)
        ax=fig.add_subplot(111)
    else:
        fig = plt.gcf()
        n_l_oplot = np.shape(ax.lines)[0]
        if n_l_oplot > 2:
            n_l_oplot=2
    fig.text(0.01, 0.01, title)

    if ax.get_xlabel()=='':
        flag_label=1

    if col[0]!='':
        colours=col
    else:
        colours_old = ['k', 'g', 'b', 'r', 'c']
        colours = ['k', 'r', 'b', 'g', 'c']
    style = styles[n_l_oplot]
    if n_lines==1:
        ax.plot(data[0], data[1], label=str(data_labels[0]), linewidth=3, color=colours[0], linestyle=style)        
    else:
        for i in range(n_lines):
            ax.plot(data[0], data[i+1], label=str(data_labels[i]), linewidth=3, color=colours[i], linestyle=style)

    if ylim!=0:
        ax.set_ylim(ylim)

    # ADJUST SUBPLOT IN FRAME
    plt.subplots_adjust(top=0.95,bottom=0.12,left=0.15,right=0.95)
    if flag_label==1:
        limit_labels(ax, xlabel, ylabel, title='')

    if data_labels[0]!='':
        ax.legend(loc='best')
    fig.tight_layout()
    plt.show()

    if fname !='':
        plt.savefig(fname, bbox_inches='tight', dpi=dpi)        

def _plot_1d(x, y=0, xlabel='', ylabel='', Id='', title='', label='',ax=0, hist=0, multip=0, fname='', color='', ylim=[0,0], ls='-'):
    """
    Hidden method for 1D plotting
    """
    common_style()

    figsize=[8,8]; flag_label=1
    flag_label=0
    if ax==0:
        # Defining figure and ax
        fig = plt.figure(figsize=figsize)
        #fig.text(0.01, 0.01, Id)
        ax  = fig.add_subplot(111)
        flag_label=1
    else:
        fig = plt.gcf()
        if multip==0:
            flag_label=0
    if hist!=0:
        if color=='':
            color='k'
        ax.hist(x,bins=30, color=color, linestyle=ls, label=label, histtype='step', lw=2.3)
    elif color=='':
        ax.plot(x,y, lw=2.3, label=label, linestyle=ls)
    else:
        ax.plot(x,y, lw=2.3, label=label, color=color, linestyle=ls)

    if ylim[0]!=0 or ylim[-1]!=0:
        ax.set_ylim(ylim)

    if flag_label==1 :
        limit_labels(ax, xlabel, ylabel, title)
        
    fig.tight_layout()
    plt.show()
    if fname !='':
        plt.savefig(fname, bbox_inches='tight', dpi=dpi)

        
def _plot_2d(x, y, xlabel='', ylabel='', dist=0, Id='', title='', wallxy=0, wallrz=0, surf=0, R0=0, ax=0, \
             scatter=0, hist=0, multip=0, xlim=0, ylim=0, fname='', cblabel='', lastpoint=1):
    """
    Hidden method to plot the 2D distribution
    wall: set to 1 if wall needed to plot (i.e. RZ function)
    """
    common_style()
    figsize=[8,6]; flag_label=1

    if ax==0:
        if wallrz!=0:
            figsize=[6,7]
        # Defining figure and ax
        fig = plt.figure(figsize=figsize)
        #fig.text(0.01, 0.01, Id)
        ax  = fig.add_subplot(111)
    else:
        ax=ax
        if multip==0:
            flag_label=0
        fig = plt.gcf()
        flag_label=0

    # Setting CB direction
    if len(fig.axes)==1:
        or_cb = 'vertical'
    else:
        or_cb = 'horizontal'

    #Doing the actual plot
    if np.mean(scatter)!=0:
        pp=ax.scatter(x, y, 40)#, c=scatter)
        #plt.colorbar(pp, ax=ax, orientation = or_cb)
    elif np.mean(dist)!=0:
        x,y = np.meshgrid(x,y)
        CS  = ax.contourf(x,y, dist, 20,  cmap=my_cmap, pad=2)
        cbar = fig.colorbar(CS, ax=ax, orientation=or_cb, shrink=1)
        cbar.ax.set_title(cblabel)     
        plt.setp(cbar.ax.get_yticklabels()[-1], visible=False) 
    elif hist != 0:
        #range = [[-3.14, 3.14],[-3.14, 3.14]]
        # ax.hist2d(x, y, bins=100, range=range, cmap=my_cmap)
        h=ax.hist2d(x, y, bins=100, cmap=my_cmap)
        cbar =fig.colorbar(h[3], ax=ax)
        cbar.ax.set_title(cblabel)
        if lastpoint==0:
            cbar.ax.set_yticklabels(cbar.ax.get_yticklabels()[0:-1])
    else:
        hb = ax.hist2d(x, y, bins=100, cmap=my_cmap)
        fig.colorbar(hb[3], ax=ax, orientation=or_cb)

    #Checks for wall and plots it	
    if wallrz != 0:
        ax.plot(wallrz[0], wallrz[1], 'k', linewidth=3)
        ax.axis('equal')
    elif wallxy != 0:
        rmin = np.min(wallxy[0])
        rmax = np.max(wallxy[0])
        circlemin = plt.Circle((0,0), rmin, color='k', fill=False, linewidth=3)
        circlemax = plt.Circle((0,0), rmax, color='k', fill=False, linewidth=3)
        ax.add_artist(circlemin); ax.add_artist(circlemax)
        ax.axis('equal')
        lims = [-rmax*1.1, rmax*1.1]
        xlim=lims; ylim=lims
    # Checks for magnetic axis XY plot
    if R0!=0:
        circle1 = plt.Circle((0, 0), R0, color='r', fill=False, linestyle='--')      
        ax.add_artist(circle1)

    #Checks for magnetic surfaces and plots them
    if surf!= 0:
        try:
            llines = [0.2, 0.4, 0.6, 0.8, 1.0]
            CS = ax.contour(surf[0], surf[1], surf[2], llines, colors='k')
            plt.clabel(CS, inline=1, fontsize=10) 
        except:
            print("Impossible to plot RZ surfaces")

    #Axes limits
    if ylim!=0:
        ax.set_ylim(ylim)
    if xlim!=0:
        ax.set_xlim(xlim)

    if flag_label == 1:
        limit_labels(ax, xlabel, ylabel, title)

    fig.tight_layout()
    plt.show()
    if fname !='':
        plt.savefig(fname, bbox_inches='tight', dpi=dpi)


def _plot_pie(x, lab, Id='', title='', ax=0, fname=''):
    """
    Hidden method to plot a pie chart
    """
    common_style()
	
    figsize=[8,8]; flag_label=1
    
    if ax==0:
        # Defining figure and ax
        fig = plt.figure(figsize=figsize)
        fig.text(0.01, 0.01, Id)
        ax  = fig.add_subplot(111)
    else:
        ax=ax
        flag_label=0
        fig = plt.gcf()

    #doing the actual plot
    plt.pie(x, labels=lab)
    
    if flag_label==1 :
        ax.axis('equal')
        
    plt.show()
    if fname!='':
        plt.savefig(fname, bbox_inches='tight', dpi=dpi)


def _plot_RZsurf(R, z, RZ, ax, surf=[0]):
    if surf[0]==0:            
    	CS = ax.contour(R, z, RZ, [0.2, 0.4, 0.6, 0.8, 1.0], colors='k')
        #plt.clabel(CS, inline=1, fontsize=10)
    else:
        CS = ax.contour(R, z, RZ, surf, colors='k')
        plt.clabel(CS, inline=True, fontsize=10, manual=True)
    return


def _cumulative_plot(x,y,labels, xlabel, ylabel, col, ax=0,  title=''):
    common_style()
    if ax==0:
        f  = plt.figure()
        ax = f.add_subplot(111)
        f.text(0.01, 0.01, title)
    else:
        ax=ax
    tmpy=np.zeros(len(x))
    for i, el in enumerate(y):
        tmpy+=el
        ax.plot(x,tmpy, col[i], lw=2.5, label=labels[i])
        ax.fill_between(x, tmpy, tmpy-el, color=col[i])
            
    limit_labels(ax, xlabel, ylabel)
    plt.tight_layout()
