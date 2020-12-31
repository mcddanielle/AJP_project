#!/usr/bin/env python3


#################################################################
'''
Supplementary Material for
Molecular dynamics simulation of synchronization in driven particles
American Journal of Physics

Danielle McDermott
Tiare Guerrero
'''
#################################################################
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#set default font of plot
plt.rc('font', size=24)

##################################################################
def landscape_force(F0,y,period):
    '''landscape potential is equation (1) in the manuscript, while force is equation (5) - for now anyway

    U(y) = U0 cos(Np*pi*y/L)

    F_y = -dU/dy = -U0/period * sin(Np*pi*y/L)

    potential minima occur at U(y) = -U0, or cos(Np*pi*y/L) = -1
    meaning Np*pi*y/L = m pi, where m is odd
    
    Required Parameters
    period: float 
    the spatial period of the washboard potential, should be the system length divided by the number of minima  (L/Np)

    y: float
    the location in the y-direction of the individual particle
    '''
    
    return -F0*np.sin(np.pi*y/period)

###################################################################
def external_drive(F_AC,F_DC,frequency,time):
    '''
    '''

    F_ext = F_DC + F_AC*np.sin(2*np.pi*frequency*time) #

    #if time % 1000:
    #    print(time,F_ext)
    
    return F_ext


####################################################################
def ramp_dc_force(F_DC,time,drop,F_DC_incr):
    '''
    '''
    F_DC_max = 0.5
    
    if time % drop == 0 and F_DC < F_DC_max:
        F_DC += F_DC_incr

    return F_DC

#####################################################################
def md_step(y, dt, time, F_DC):
    '''
    '''

    F_AC = 0.05
    #F_DC = 0.005 #5 #increment slowly
    period = 1.825
    F0 = 0.1 #np.pi/period
    frequency = 0.01
    SX = 36.5
    SY = 36.5
    
    #integrate the equation of motion \eta v = F^net (eta=1)

    #reset vy for every timestep!  
    vy = landscape_force(F0,y,period)

    #driving force
    vy += external_drive(F_AC,F_DC,frequency,time)

    #calculate the new position
    y += vy*dt

    #check periodic boundary conditions
    if y > SY:
        y -= SY
    elif y < 0:
        y += SY
    
    return y, vy


#-------------------------------------------------------------------
if __name__ == "__main__":

    #define the constants (parameters and initial conditions)
    dt = 0.001 #
    Sy=36.5
    vy = 0
    time = 0.0
    period = 1.825

    #start the particle in a local minima
    y = period #Sy/2

    F_DC = 0
    F_DC_incr = 0.01
    F_DC_max = 0.1
    drop = 4000
    
    #integer time steps
    maxtime=300000
    writemovietime=1000
    
    #define arrays to hold data as a function of time
    array_length = int(maxtime/writemovietime)
    y_data = np.zeros(array_length)
    FDC_data = np.zeros(array_length)
    vy_data = np.zeros(array_length)
    time_data = np.zeros(array_length)

    #loop through the integer time steps in the simulation
    for int_time in range(0,maxtime):

        if int_time > 0 and F_DC < F_DC_max: 
            F_DC = ramp_dc_force(F_DC,int_time,drop,F_DC_incr)
            #print(F_DC)

        #apply the force calculations for the current position/time
        #note vy is not necessary 
        y, vy = md_step(y, dt, time, F_DC)

        if int_time % writemovietime == 0:
            #print(int_time)
            #update the data
            i = int(int_time/writemovietime)
            y_data[i] = y
            vy_data[i] = vy
            FDC_data[i] = F_DC
            time_data[i] = time

        #update the time (simulation units)
        time += dt

    #calculate the driving force from F_DC and time
    F_AC = 0.05
    freq = 0.01
    F_drive = FDC_data + F_AC*np.sin(2*np.pi*freq*time_data)

    #max occurs at
    # 2*np.pi*freq*time_data = m*pi/2, where m = 1, 5, 9, 13, ...
    # time = m/(4*freq)

    fig = plt.figure(figsize=(8,8))
    gs=gridspec.GridSpec(2,1)
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
    ax2 = fig.add_subplot(gs[1,0]) #,sharex=ax1)  #scatter plot of particles
    
    #plot the data
    ax1.plot(time_data,FDC_data,label="F$_{DC}$")
    ax1.plot(time_data,F_drive,label="F$^D$(t)",lw=5)
    ax1.legend(loc=4,fontsize=20,borderaxespad=0.1,frameon=0) #,labelspacing=0.2
    ax2.plot(time_data,y_data/period,lw=5) #,'o--') #,markevery=100)
    #plt.xlim(155,157)
    ax2.set_xlabel("time")
    ax2.set_ylabel(r"y/$\lambda$")
    ax1.set_ylabel(r"F$^D$(t)")

    ax1.set_xticks([])

    ax1.set_xlim(0,time_data[-1]+1)
    ax2.set_xlim(0,time_data[-1]+1)

    #add horizontal lines for potential minima
    for i in [1,3,5,7]:
        #ax2.hline(i)
        #, label = "U(y) = -U$_0$")
        ax2.axhline(y=i, color = "black", linestyle = "--")

    for i in [1,5,9,13]:
        #peaks of the FD curve
        ax2.axvline(x=i/(4*freq), color = "black", linestyle = "--") 
        ax1.axvline(x=i/(4*freq), color = "black", linestyle = "--") 
     
    plt.tight_layout(pad=0.1)
    plt.savefig("single_particle.pdf")
