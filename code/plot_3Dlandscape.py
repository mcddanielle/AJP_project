#!/usr/bin/env python3

#numeric and basic plotting libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from matplotlib import cm
import matplotlib.cm as cm

import matplotlib.gridspec as gridspec

from mpl_toolkits.mplot3d import Axes3D
#import mpl_toolkits.mplot3d.Axes3D as Axes3D


#set default font of plot
plt.rc('font', size=20)

##########################################################
#ADD CONTOUR PLOT
##########################################################
def add_contour(ax,L,N,corrugated = True):
    '''
    Hardwired to color in the quasi1D potential to contain 
    the particles in a trough.  
    Can also add the washboard/corrugated substrate.

    Required Arguments

    Optional Arguments:

    corrugated (default = True)  
    Adds the washboard in the y-direction.  
    Hardwired for a single parameter set.    
    '''

    #set up plot of 1 row, 3 columns to make the first plot smaller
    gs=gridspec.GridSpec(1,3)
    
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
    ax2 = fig.add_subplot(gs[0,1:], projection='3d')  #scatter plot of particles

    a_p = L/N

    #assuming Tiare's trough system, so we won't want to cover the entire range
    X = np.arange(0, L/2, 0.1)
    Y = np.arange(0, L, 0.1)
    X, Y = np.meshgrid(X, Y)

    Z_mag = 0.1 # set by what "looks good"

    #create a narrow confining channel for 2D systems
    if 0:
        Z = Z_mag*np.cos(2*np.pi*X/L)
    
        #create washboard
    elif corrugated == True:
        Z = Z_mag*np.cos(2*np.pi*(Y)/a_p) 

    cmap=cm.coolwarm #_r

    #alphs is the degree of transparency, again, set by what looks good.

    #if 0:
    cset = ax1.contourf(X, Y, Z, cmap=cmap,alpha=0.5)
    #else:
    cset = ax2.plot_surface(X,Y,Z,rstride=1,cstride=1,linewidth=0,cmap=cmap,alpha=0.5)

    #setboxaspect doesn't work in my version of matplotlib
    #ax2.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))  # aspect ratio is 1:1:1 in data space

    #the following doesn't seem to work at all
    ax2.auto_scale_xyz([0,L/2], [0, L], [-Z_mag, Z_mag])

    
    ax2.set_zticks([])
    #if stripes:

    ax2.w_xaxis.set_pane_color((1.0,1.0,1.0,1.0))
    ax2.w_yaxis.set_pane_color((1.0,1.0,1.0,1.0))
    ax2.w_zaxis.set_pane_color((1.0,1.0,1.0,1.0))

    ax2.set_xlim(-1, L/2)
    ax2.set_ylim(-1, L)
    ax2.set_zlim(-2*Z_mag, 2*Z_mag)

    ax1.set_xlim(0, L/2)
    ax1.set_ylim(0, L)
    #ax.view_init(35, 60)
    #ax.dist = 12

    #ax1.set_xlim(15,20)
    #ax1.set_ylim(15,20)

    for ax in [ax1,ax2]:
        ax.set_xlabel(r"X")
        ax.set_ylabel(r"Y",rotation='horizontal',ha='right')

        ax.set_xticks([0,L/2])
        ax.set_yticks([0,L])

        ax.set_xticklabels(['0','L/2'])
        ax.set_yticklabels(['0','L'])

    #ax1.set_xticks([])
    #ax1.set_yticks([])
    return ax1, ax2

################################################################
################################################################
################################################################

if __name__ == "__main__":

    #system size
    Sx=36.5
    Sy=36.5


    fig = plt.figure(figsize=(12,6))
    '''
    if 1:
        #ax1 = fig.add_subplot(121)
        #ax2 = fig.add_subplot(122, projection='3d')
        ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
        ax2 = fig.add_subplot(gs[0,1:], projection='3d')  #scatter plot of particles
    else:
        ax1 = Axes3D(fig)
    '''
    
    #set up landscape
    n_corr = 4  #number of periods
    corr = True
    ax1, ax2 = add_contour(fig,Sy,n_corr,corrugated = corr)

    #fig2 = plt.figure()
    #ax2 = Axes3D(fig2)

    #place single particle on plot
    xp=Sx*0.25
    yp=Sy*0.5
    scatter1=ax1.scatter(xp,yp,edgecolor='k',s=50)

    #add arrow for driving force

    #annotate variables on plot

    #configure and save the image
    fig.tight_layout() #h_pad=-0.5,w_pad=1.0,pad=0.5)
    image_test_name = "landscape.png"
    plt.savefig(image_test_name)

    #sys.exit()
