#!/usr/bin/env python3

#numeric and basic plotting libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from matplotlib import cm
import matplotlib.cm as cm

import matplotlib.gridspec as gridspec


#set default font of plot
plt.rc('font', size=24)

##########################################################
#ADD CONTOUR PLOT
##########################################################
def add_contour(ax,ax2,L,N): 
    '''
    Hardwired to color in the quasi1D potential to contain 
    the particles in a trough.  
    Can also add the washboard/corrugated substrate.

    Required Arguments

    Optional Arguments:

    Adds the washboard in the y-direction.  
    Hardwired for a single parameter set.    
    '''

    a_p = L/N

    #assuming Tiare's trough system, so we won't want to cover the entire range
    X = np.arange(0, L/2, 0.1)
    Y = np.arange(0, L, 0.1)
    X, Y = np.meshgrid(X, Y)

    Z_mag = 0.1 # set by what "looks good"

    #create a narrow confining channel for multiparticle 2D systems
    if 0:
        Z = Z_mag*np.cos(2*np.pi*X/L)
    
        #create washboard

    Z = Z_mag*np.cos(2*np.pi*(Y)/a_p) 

    #blue are minima, red are maxima
    cmap=cm.coolwarm #_r

    #alpha is the degree of transparency, again, set by what looks good.

    cset = ax.contourf(X, Y, Z, cmap=cmap,alpha=0.5)

    ax2.plot(Z,Y)
    ax2.set_xlabel(r"U(y)",labelpad=-20)
    #ax2.set_ylabel(r"y",rotation='horizontal',ha='right')

    ax.set_xlim(0, L/2)
    ax.set_ylim(0, L)

    ax.set_xlabel(r"x",labelpad=-20)
    ax.set_ylabel(r"y",rotation='horizontal',ha='right',labelpad=-15)

    ax.set_xticks([0,L/2])
    ax.set_yticks([0,L])

    ax.set_xticklabels(['0','L/2'])
    ax.set_yticklabels(['0','L'])

    #ax2.set_xlim(-Z_mag*1.2,Z_mag*1.2)
    ax2.set_xticks([-Z_mag,Z_mag])
    ax2.set_xticklabels([r'-U$_0$','U$_0$'])
    #ax2.set_yticks([])
    return 

################################################################
################################################################
################################################################

if __name__ == "__main__":

    #system size
    Sx=36.5
    Sy=36.5


    fig = plt.figure(figsize=(8,6))

    #set up plot of 1 row, 3 columns to make the first plot smaller
    gs=gridspec.GridSpec(1,2)
    
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
    ax2 = fig.add_subplot(gs[0,1],sharey=ax1)  #scatter plot of particles
    
    #set up and plot landscape
    n_corr = 3  #number of periods

    add_contour(ax1,ax2,Sy,n_corr) #,corrugated = corr)

    #place single particle on plot
    xp=Sx*0.25
    yp=Sy*0.41
    Z_mag=0.1
    zp=Z_mag*np.cos(2*np.pi*yp/(Sy/n_corr))

    #slope at zp,yp:
    slope = -Z_mag*2*np.pi/(Sy/n_corr)*np.sin(2*np.pi*yp/(Sy/n_corr))

    #slope = dy/dz, dz=1, zf-zp=Z_mag, yf-yp=slope
    zf=zp-Z_mag/1.8+0.02
    yf=yp-slope*Z_mag/1.8+0.6
    

    scatter1=ax1.scatter(xp,yp,edgecolor='k',s=200)
    scatter2=ax2.scatter(zp,yp,edgecolor='k',s=200)

    ax2.annotate(r" $\vec{F}^{landscape}$", xytext=(zp, yp+0.5), xy=(zf, yf), arrowprops=dict(facecolor='black', shrink=0.05)) #arrowstyle="<-",

    ax1.annotate(r" $\vec{F}^{landscape}$", xytext=(xp, yp+4), xy=(xp, yp), arrowprops=dict(facecolor='black', arrowstyle="<-"),ha="center")

    ax1.text(0.03,0.94,"(a)",transform = ax1.transAxes)
    ax2.text(0.03,0.94,"(b)",transform = ax2.transAxes)
    
    #plot the "side view" of the potential

    #add arrow for driving force

    #annotate variables on plot

    #configure and save the image
    fig.tight_layout(pad=0.5) #h_pad=-0.5,w_pad=1.0,pad=0.5)
    image_test_name = "landscape.png"
    plt.savefig(image_test_name)

    #sys.exit()
