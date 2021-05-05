#!/usr/bin/env python3

#################################################################
'''
Supplementary Material for
Molecular dynamics simulation of synchronization in driven particles
American Journal of Physics

Solutions to Problems 2,3,5,6

Danielle McDermott
Tiare Guerrero
'''
#################################################################
import numpy as np
import math, sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

from scipy.optimize import curve_fit

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from mpl_toolkits.axes_grid.inset_locator import (inset_axes,
                                                  InsetPosition,
                                                  mark_inset)


#set default font of plot
plt.rc('font', size=24)

#main MD routines
import MD_colloid

#############################################################
#function that curve_fit will call/fit for problem 6
#############################################################
def func(x, a, b, c):
    '''
    function to use in scipy optimize curve_fit call

    required arguments
    x is an independent variable
    a,b,c are fit parameters which are returned by value
    
    returns:
    popt : array, Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized
    pcov : 2d array, The estimated covariance of popt.
    '''
    
    return a*((x+c)**(-b))
#-------------------------------------------------------------------

if __name__ == "__main__":

    #will pass the parameters dict to subroutines to read constants
    parameters = MD_colloid.set_parameters()

    #select which figure in the AJP you would like to make
    #figures 2a, 2b, 2c, 5, 6
    solve_problem = "6"

    parameters['filename']="problem_%s.pdf"%(solve_problem)

    letter=["(a)","(b)","(c)","(d)","(e)","(f)","(g)"]
        
    #--------------------------------------------------------------
    #modify either frequency (2b) or F^ac amplitude (2a) or F^ac = 0
    #--------------------------------------------------------------
    if solve_problem == "2a" or solve_problem == "2b":
        
        parameters['dt'] = 0.1 #timestep in simulation units
        parameters['writemovietime']=1   #interval to write data to arrays for plotting

        if solve_problem == "2a":
            parameters['maxtime']=3000        #total time steps in simulation

            #we will sweep through the following independent variable
            F_AC = [0.2, 0.3, 0.4]
            ind_var = F_AC
            str_iv=" $F_{\mathrm{ac}}$=%1.2f"
            #parameters['freq'] = 0.1
            #parameters['F_DC'] = 0.1
            
        else:
            parameters['maxtime']=4500        #total time steps in simulation

            #we will sweep through the following independent variable

            frequency = [0.1,0.015,0.01, 0.005, 0.002, 0.001]
            ind_var = frequency
            str_iv=" f=%1.3f"

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
            if solve_problem == "2b":
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
    #make shapiro steps with a range of F^{dc}
    #this is the place to include the strictly DC force...
    #--------------------------------------------------------------
    if solve_problem == "2c": 

        #sweep through a variety of dc values

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

        ax1.set_xlim(0,Fdc_data[-1])
        ax1.set_ylim(avg_vy_data[0]-0.1,avg_vy_data[-1])

        #plt.xlim(155,157)
        ax1.set_xlabel("F$_{\mathrm{dc}}$")
        ax1.set_ylabel(r"$\langle$v$_y \rangle / \lambda f $")
        ax1.legend(loc='best',fontsize=20,numpoints=3,borderpad=0.1,handletextpad=0.2,title="F$_{\mathrm{ac}}$")
                    
        plt.tight_layout(pad=0.2)
        plt.savefig(parameters['filename'])

    #--------------------------------------------------------------
    #calculate the approximate Reynolds number
    #at standard temperature and pressure, converted to SI units (i.e. m-kg-sec)

    if solve_problem == "3":

        D = 1e-6 #micrometer
        rho = 0.9982 * 1e-3 * 1e6 #g/cm^3  1g=10^-3kg, 1/cm^3*(100cm/1m)^3
        v = 1e-6 #micrometer / sec

        #"kinematic viscosity" - not the right units
        #eta = 0.01 #cm^2 / sec

        #"dynamics viscosity" - same units as Taylor
        eta = 1.0016e-3 #mPa *sec = mN/m^2 *sec


        #unit check:  m*m/sec*kg/(m^3)*1/(N/m^2*sec)
        #blech! 
        #unit check:  kg/m * m^2 / (N*sec^2) 
        #unit check:  m*kg / (N*sec^2) = m*kg / (kg*m s^2 /s^2)
        #kg cancels!  m cancels! s^2 cancels!
        #unitless Reynold's number


        R = D*v*rho/eta
        print("The Reynolds number is: ",R)
        
    #--------------------------------------------------------------
    #study Brownian motion for increasing temperature
    #--------------------------------------------------------------
    if solve_problem == "5":

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
        for temp in [4.0, 5.0, 6.0]:

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
        
    #--------------------------------------------------------------
    #Problem 6
    #Fit depinning exponent using data from Fig. 3
    #--------------------------------------------------------------
    if solve_problem == "6":

        #substrate height
        parameters['AP'] = 0.1

        #sweep through a variety of dc values
        delta_Fdc = 0.001
        Fdc_max=0.2+delta_Fdc

        #DC signal
        parameters['F_AC'] = 0.0

        #set up arrays
        max_value = int(np.ceil(Fdc_max/delta_Fdc))
        avg_vy_data = np.zeros(max_value)
        Fdc_data = np.arange(0,Fdc_max,delta_Fdc)

        #make the figure 
        fig = plt.figure(figsize=(7,4.5))
        gs=gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles

        #use m to find the integer range of fit values
        #above the critical force
        m=0
        
        #run the MD sim for each value of FD, this takes a while...
        for i in range(len(avg_vy_data)):

            #sweep through a range of dc values
            parameters['F_DC'] = Fdc_data[i] 
                
            #run a single MD simulation for a set of parameters
            avg_vy_data[i] = MD_colloid.single_particle(parameters,plot="y-velocity")

            #identify depinning to select the data to fit
            #where the vy is non-zero
            if m == 0 and (avg_vy_data[i])/(parameters['freq']*parameters['period']) > 0.18:
                #step back to the last step before depinning
                m = i-1
                print("fitting arrays for Fdc > ",Fdc_data[i])
                
        #normalize <vy> so you can count step number
        avg_vy_data /= (parameters['freq']*parameters['period'])
            
        #plot <vy> vs F_dc
        ax1.plot(Fdc_data,avg_vy_data,"o",linewidth=2,label="MD data")

        #curve fit
        popt, pcov = curve_fit(func, Fdc_data[m:-1], avg_vy_data[m:-1])
        perr = np.sqrt(np.diag(pcov))

        #plot curve fit
        ax1.plot(Fdc_data[m:-1], func(Fdc_data[m:-1], *popt), 'r-',
                 label=r'fit: $\beta$=%5.3f$\pm$%5.3f, F$_c$=%6.4f$\pm$%6.4f' % (-popt[1],perr[1],-popt[2],perr[2]))

        #make a nice plot
        ax1.set_xlim(0,Fdc_data[-1])
        ax1.set_ylim(avg_vy_data[0]-0.1,avg_vy_data[-1])
        ax1.set_xlabel("F$_{\mathrm{dc}}$")
        ax1.set_ylabel(r"$\langle$v$_y \rangle / \lambda f $")
        ax1.legend(loc=2,fontsize=16,numpoints=3,borderpad=0.1,handletextpad=0.2)                    
        plt.tight_layout(pad=0.2)
        plt.savefig(parameters['filename'])
        
    sys.exit()


