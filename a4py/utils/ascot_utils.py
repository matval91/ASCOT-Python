import matplotlib.pyplot as plt
from matplotlib import ticker, colors
import numpy as np

def define_colors():
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

    return colours, colours_old, styles, my_cmap, dpi


colours, colours_old, styles, my_cmap, dpi= define_colors()


def common_style():
    """plotting style

    Defining common, nice plotting styles using plt.rc
    
    Parameters:
        None
    Arguments:
        None
    """      
    plt.rc('font', weight='bold')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('axes', labelsize=30, labelweight='normal', titlesize=24)
    plt.rc('figure', facecolor='white')
    plt.rc('legend', fontsize=20)

def limit_labels(ax, xlabel='', ylabel='', title='', M=5):
    """plotting style

    Limiting labels of the axis to 4 elements and setting grid
    
    Parameters:
        | ax (axis object) : axis object where to limit labels
        | xlabel (str) : xlabel of the plot (default '')
        | ylabel (str) : ylabel of the plot (default '')
        | title (str)  : title of the plot (default '')
        | M (int)      : number of major ticks in x and y axis (default 5)
    Arguments:
        None
    """   

    #==============================================
    # SET TICK LOCATION
    #==============================================
    
    # Create your ticker object with M ticks
    yticks = ticker.MaxNLocator(M)
    xticks = ticker.MaxNLocator(M)
    ax.yaxis.set_major_locator(yticks)
    ax.xaxis.set_major_locator(xticks)
    #==============================================
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid('on', alpha=0.6)
    
    #Removing first point of y-axis
    plt.setp(ax.get_yticklabels()[0], visible=False) 


def plot_article(n_lines, data, data_labels, xlabel='', ylabel='',\
                 title='', ax=0, ylim=0, fname='', col=['']):
    """
        calls plot_multilines (done this way in order not to redefine everything)
    """
    plot_multilines(n_lines, data, data_labels, xlabel, ylabel, title, ax, ylim, fname, col)
    
    
def plot_multilines(n_lines, data, data_labels, xlabel='', ylabel='', \
                    title='', ax=0, ylim=0, fname='', col=['']):
    """plotting multiple lines

    Plot multiple lines together on same plot
    
    Parameters:
        | n_lines (int)         : number of lines to plot
        | data (array)          : data to plot. shape=[n_lines,xpoints]. usually data[0]
                                is x axis
        | data_labels (array)   : array with labels of each line
        | xlabel (str)          : xlabel of the plot (default '')
        | ylabel (str)          : ylabel of the plot (default '')
        | title (str)           : title of the plot (default '')
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)
        | ylim (array)          : set y limit (default=0, no limit set) 
        | fname (str)           : file name where to save (default '')
        | col (array)           : array with the colors, if wanted different from default
    Arguments:
        None
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
        colours = ['k', 'r', 'b', 'g', 'c']
    style = styles[n_l_oplot]
    style='-'
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

def _plot_1d(x, y=0, xlabel='', ylabel='', title='', label='',ax=0, hist=0, \
             ylim=[0,0], fname='', color='', ls='-'):
    """plotting 1D

    Private method for 1D plotting (x vs y)
    
    Parameters:
        | x (array)             : independent variable
        | y (array)             : dependent variable (default=0, i.e. plots hist)
        | xlabel (str)          : xlabel of the plot (default '')
        | ylabel (str)          : ylabel of the plot (default '')
        | title (str)           : title of the plot (default '')
        | label (array)         : label of the line
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)
        | hist (int)            : if !=0, plots the histogram (default=0)
        | ylim (array)          : set y limit (default=0, no limit set) 
        | fname (str)           : file name where to save (default '')
        | color (str)           : color of line, if wanted different from default
        | ls (str)              : linestyle (default='-')
    Arguments:
        None
    """ 
    common_style()

    figsize=[8,8]
    if ax==0:
        # Defining figure and ax
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111)
    else:
        fig = plt.gcf()
    
    if color=='': color='k';
        
    if hist!=0:
        ax.hist(x,bins=30, color=color, linestyle=ls, label=label, histtype='step', lw=2.3)
    ax.plot(x,y, lw=2.3, label=label, color=color, linestyle=ls)

    if ylim[0]!=0 or ylim[-1]!=0:
        ax.set_ylim(ylim)
      
    fig.tight_layout()
    plt.show()
    if fname !='':
        plt.savefig(fname, bbox_inches='tight', dpi=dpi)

        
def _plot_2d(x, y, xlabel='', ylabel='', dist=0, title='', wallxy=0, wallrz=0, surf=0, R0=0, ax=0, \
             scatter=0, hist=0, xlim=0, ylim=0, fname='', cblabel='', lastpoint=1):
    """plotting 2D

    Private method for 2D plotting (e.g. distribution functions!)
    
    Parameters:
        | x (array)             : independent variable
        | y (array)             : independent variable 
        | xlabel (str)          : xlabel of the plot (default '')
        | ylabel (str)          : ylabel of the plot (default '')
        | dist (array)          : if plotting a distribution, give the values here.
                                    (x,y) will be the meshgrid (default 0)
        | title (str)           : title of the plot (default '')
        | wallxy (array)        : 2D array with [R,Z] of the wall. Will plot circles with
                                    min(R) and max(R) (default 0)
        | wallrz (array)        : 2D array with [R,Z] of the wall. Will plot RZ wall section
                                    (default 0)
        | surf (array)          : 2D array with [R,Z, surf] of the flux surfaces. 
                                    Will plot flux surfaces wrt R,Z (default 0)      
        | R0 (float)            : plots circle with value R0 (default 0)        
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)
        | scatter (int)         : if !=0, plots scatter of (x,y) (default=0)        
        | hist (int)            : if !=0, plots the histogram (default=0)
        | xlim (array)          : set x limit (default=0, no limit set) 
        | ylim (array)          : set y limit (default=0, no limit set) 
        | fname (str)           : file name where to save (default '')
        | cblabel (str)         : label of colorbar (default='')
        | lastpoint (str)       : if 0, doesn't plot last point in colorbar (default=1)
    Arguments:
        None
    """     
    common_style()
    figsize=[8,6]; flag_label=1

    if ax==0:
        if wallrz!=0:
            figsize=[6,7]
        # Defining figure and ax
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111)
    else:
        ax=ax
        fig = plt.gcf()
        flag_label=0

    # Setting CB direction
    or_cb = 'horizontal'
    if len(fig.axes)==1:
        or_cb = 'vertical'
    #Doing the actual plot
    if np.mean(scatter)!=0:
        ax.scatter(x, y, 40)#, c=scatter)
    elif np.mean(dist)!=0:
        x,y = np.meshgrid(x,y)
        CS  = ax.contourf(x,y, dist, 20,  cmap=my_cmap, pad=2)
        cbar = fig.colorbar(CS, ax=ax, orientation=or_cb, shrink=1)
        cbar.ax.set_title(cblabel)     
        plt.setp(cbar.ax.get_yticklabels()[-1], visible=False) 
    elif hist != 0:
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


def _plot_pie(x, lab, title='', ax=0, fname=''):
    """plotting pie

    Private method for pie plotting
    
    Parameters:
        | x (array)             : variables of pie plot
        | lab (array)           : labels of plot
        | title (str)           : title of the plot (default '')      
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)
        | fname (str)           : file name where to save (default '')

    Arguments:
        None
    """   
    common_style()
	
    figsize=[8,8]; flag_label=1
    
    if ax==0:
        # Defining figure and ax
        fig = plt.figure(figsize=figsize)
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
    """plotting flux surf

    Private method for plotting surf
    
    Parameters:
        | R (array)             : major radius for plot
        | z (array)             : z axis for plot
        | RZ (str)              : 2D flux surf grid      
        | ax (axis object)      : axis object where to plot
        | surf (arr)            : define which surfaces to plot (default [0])

    Arguments:
        None
    """  
    if surf[0]==0:            
        	CS = ax.contour(R, z, RZ, [0.2, 0.4, 0.6, 0.8, 1.0], colors='k')
    else:
        CS = ax.contour(R, z, RZ, surf, colors='k')
        plt.clabel(CS, inline=True, fontsize=10, manual=True)
    return


def _cumulative_plot(x,y,labels, xlabel, ylabel, col, ax=0,  title=''):
    """plotting cumulative

    Private method for plotting cumulative lines (e.g. power balance with transp)
    
    Parameters:
        | x (array)             : independent variable
        | y (array)             : dependent variables 2D array with [n, xpoints]
        | labels (str)          : labels of the different lines 
        | xlabel (str)          : xaxis label
        | ylabel (arr)          : yaxis label
        | col (str)             : color of lines
        | ax (axis object)      : axis object where to plot (default=0, creates a new figure)
        | title (str)           : title of the plot (default '')
    Arguments:
        None
    """  
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
