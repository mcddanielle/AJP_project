#!/usr/bin/env python3

#################################################################
'''
Supplementary Material for
Molecular dynamics simulation of synchronization in driven particles
American Journal of Physics

Makes Figures 2-9 using a flag to select

Danielle McDermott
Tiare Guerrero
'''
#################################################################
import numpy as np
import math, sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from mpl_toolkits.axes_grid.inset_locator import (inset_axes,
                                                  InsetPosition,
                                                  mark_inset)


#set default font of plot
plt.rc('font', size=24)

#main MD routines
import MD_colloid


#-------------------------------------------------------------------
if __name__ == "__main__":

    #will pass the parameters dict to subroutines to read constants
    parameters = MD_colloid.set_parameters()

    #select which figure in the AJP you would like to make
    #figures 2a, 2b, 2c, 5
    solve_exercise = "5"

    parameters['filename']="exercise_%s.pdf"%(solve_exercise)

    letter=["(a)","(b)","(c)","(d)","(e)","(f)","(g)"]

    #--------------------------------------------------------------
    #make shapiro steps with a range of F^{dc}
    #this is the place to include the strictly DC force...
    #--------------------------------------------------------------
    if solve_exercise == "3c": 

        #sweep through a variety of dc values
        delta_Fdc = 0.001

        #IN EXPLORING THE HIGH FREQUENCY REGIME,
        #I'M SEEING STEP HEIGHTS LIKE JUNIPER (but still fractional)  
        parameters['freq'] = 0.1
        delta_Fdc = 0.01
        F_AC = [0.1,0.2,0.3,0.4]
        Fdc_max=0.6+delta_Fdc

        max_value = int(np.ceil(Fdc_max/delta_Fdc))
        avg_vy_data = np.zeros(max_value)
        Fdc_data = np.arange(0,Fdc_max,delta_Fdc)

        #make the figure 
        fig = plt.figure(figsize=(7,4.5))
        gs=gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles

        for F in F_AC:
            parameters['F_AC'] = F

            #run the MD sim for each value of FD, this takes a while...
            for i in range(len(avg_vy_data)):

                #sweep through a range of dc values
                parameters['F_DC'] = Fdc_data[i] 
                
                #run a single MD simulation for a set of parameters
                avg_vy_data[i] = MD_colloid.single_particle(parameters,plot="y-velocity")
                
            #normalize <vy> so you can count step number
            avg_vy_data /= (parameters['freq']*parameters['period'])
            
            #plot <vy> vs F_dc
            if F == 0.0:
                ax1.plot(Fdc_data,avg_vy_data,"--",linewidth=2,label=r"%1.2f"%(F))
            else:
                ax1.plot(Fdc_data,avg_vy_data,linewidth=2,label=r"%1.2f"%(F))

            #make the phase plot locations
            if make_fig == 3 and F > 0:
                for i in [40,70,100,125]:
                    ax1.scatter(Fdc_data[i],avg_vy_data[i],c='k',zorder=10)
            
        ax1.set_xlim(0,Fdc_data[-1])
        ax1.set_ylim(avg_vy_data[0]-0.1,avg_vy_data[-1])

        #plt.xlim(155,157)
        ax1.set_xlabel("F$_{\mathrm{dc}}$")
        ax1.set_ylabel(r"$\langle$v$_y \rangle / \lambda f $")
        ax1.legend(loc='best',fontsize=20,numpoints=3,borderpad=0.1,handletextpad=0.2,title="F$_{\mathrm{ac}}$")
                    
        plt.tight_layout(pad=0.2)
        plt.savefig(parameters['filename'])
        
    #--------------------------------------------------------------
    #modify either frequency (fig5) or F^ac amplitude (fig6) or F^ac = 0
    #--------------------------------------------------------------
    if solve_exercise == "2a" or solve_exercise == "2b":
        
        parameters['dt'] = 0.1 #timestep in simulation units
        parameters['writemovietime']=1   #interval to write data to arrays for plotting


        if solve_exercise == "2a":
            parameters['maxtime']=4500        #total time steps in simulation

            #we will sweep through the following independent variable

            frequency = [0.1,0.015,0.01, 0.005, 0.002, 0.001]
            ind_var = frequency
            str_iv=" f=%1.3f"

        else:
            parameters['maxtime']=3000        #total time steps in simulation

            #we will sweep through the following independent variable
            F_AC = [0.2, 0.4]
            ind_var = F_AC
            str_iv=" $F_{\mathrm{ac}}$=%1.2f"
            #parameters['freq'] = 0.1
            #parameters['F_DC'] = 0.1

        #make the figure
        n_plots = len(ind_var)
        fig = plt.figure(figsize=(10,3*n_plots))
        gs=gridspec.GridSpec(n_plots,1)
        
        #run a single MD simulation for a set of parameters in independent variable
        for i in range(n_plots):

            #make the subplot
            ax = fig.add_subplot(gs[i,0])
            ax.text(0.7,0.08,letter[i]+str_iv%(ind_var[i]),
                    transform = ax.transAxes,backgroundcolor="white")
            parameters['axis'] = ax

            #select the independent variable
            if make_fig == 5:
                parameters['freq'] = frequency[i]
            else:
                parameters['F_AC'] = F_AC[i]

                #horizontal grid so we can count hops
                #ax.set_ylim(-0.9,20.9)
                ax.yaxis.set_minor_locator(AutoMinorLocator(1))
                ax.grid(axis='y',which='both')

            #run the MD simulation    
            MD_colloid.single_particle(parameters,plot="just_position")

            if i < n_plots-1:
                ax.tick_params(labelbottom=False)
                ax.set_xlabel("")

        #save the figure
        plt.tight_layout(pad=0.1,h_pad=-0.05)
        plt.savefig(parameters['filename'])

    #--------------------------------------------------------------
    #study Brownian motion for increasing temperature
    #--------------------------------------------------------------
    if solve_exercise == "5":

        #place particle in center to reduce hops across PBC
        parameters['y0'] = parameters['Sy']/2

        #simulating for a long time, 
        parameters['dt'] = 0.1 #timestep in simulation units

        parameters['F_DC'] = 0.0 #don't drive
        parameters['F_AC'] = 0.0 #don't drive

        #run for a very long time to see the infrequent jumps
        parameters['maxtime']=300000        #total time

        #that's a lot of data!
        parameters['writemovietime']=1000   #interval to write data 


        #make the figure
        fig = plt.figure(figsize=(10,12))
        gs=gridspec.GridSpec(3,1)
        i=0

        #sweep through temperatures
        for temp in [3.0, 3.5, 4.0]:

            #create and label the subplot
            ax = fig.add_subplot(gs[i,0])
            ax.text(0.02,0.05,letter[i],
                    transform = ax.transAxes,backgroundcolor="white")

            #Brownian motion factor as a ratio of the Ap value
            parameters['temperature'] = temp*parameters['AP']

            #run a single MD simulation for a set of parameters,
            #will be plotted on subplot
            MD_colloid.brownian_particle(parameters,ax)

            i+=1

        #save fig 7
        plt.tight_layout(pad=0.1,h_pad=-0.3)
        plt.savefig(parameters['filename'])
        
    #--------------------------------------------------------------------------

    sys.exit()


