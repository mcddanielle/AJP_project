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
import math, sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#set default font of plot
plt.rc('font', size=24)

##################################################################
def landscape_force(y,parameters):
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

    period = parameters['period']
    F0 = parameters['F0']
    
    return -F0*np.sin(np.pi*y/period)

###################################################################
def external_drive(F_DC,time, parameters):
    '''
    '''

    F_AC = parameters['F_AC']
    frequency = parameters['freq']
    
    F_ext = F_DC + F_AC*np.sin(2*np.pi*frequency*time) #

    #if time % 1000:
    #    print(time,F_ext)
    
    return F_ext


####################################################################
def ramp_dc_force(F_DC, time, parameters):
    '''
    '''
    F_DC_max = parameters['F_DC_max']
    drop = parameters['drop']
    F_DC_incr = parameters['F_DC_incr']
    
    if time % drop == 0 and F_DC < F_DC_max:
        F_DC += F_DC_incr

    return F_DC

#####################################################################
def average_velocity(time, vy, avg_vy, parameters):
    '''calculate the average velocity over a range of time values for N particles.  
    Use all vy values, not just those saved during "writemovietime.  So far just "
    written for the single particle case, where the code is not so interesting.  Subroutine 
    generally sums all vxi and vyi values across n particles and averages
    '''

    decifactor = parameters['decifactor']

    avg_vy += vy
            
    if time == decifactor:

        avg_vy /= decifactor

    return avg_vy


#####################################################################
def md_step(y, time, F_DC, avg_vy, parameters):
    '''
    '''

    dt = parameters['dt']
    SY = parameters['Sy']
    
    #integrate the equation of motion \eta v = F^net (eta=1)

    #reset vy for every timestep since ay = 0 
    vy = landscape_force(y, parameters)

    #driving force
    vy += external_drive(F_DC, time, parameters)

    #calculate the average velocity over all vy values
    avg_vy = average_velocity(time, vy, avg_vy, parameters)
    
    #calculate the new position
    y += vy*dt

    #check periodic boundary conditions
    if y > SY:
        y -= SY
    elif y < 0:
        y += SY
    
    return y, vy, avg_vy


#####################################################################
def plot_force_position_vs_time(time_data,FDC_data,y_data,p):
    '''Make Fig. 2 in AJP
    '''

    #for plot
    #calculate the driving force from F_DC and time
    F_AC = p['F_AC']
    freq = p['freq']
    period = p['period']
    
    fig = plt.figure(figsize=(8,8))
    gs=gridspec.GridSpec(2,1)
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
    ax2 = fig.add_subplot(gs[1,0]) #,sharex=ax1)  #scatter plot of particles

    #calculate from np.array rather than incrementally from subroutines
    F_drive = FDC_data + F_AC*np.sin(2*np.pi*freq*time_data)

    #max occurs at
    # 2*np.pi*freq*time_data = m*pi/2, where m = 1, 5, 9, 13, ...
    # time = m/(4*freq)

        #plot the data
    ax1.plot(time_data,FDC_data,label="F$^{dc}$")
    ax1.plot(time_data,F_drive,label="F$^d$(t)",lw=5)
    ax1.legend(loc=4,fontsize=20,borderaxespad=0.0,frameon=0,handlelength=1.5) #,labelspacing=0.2
    ax2.plot(time_data,y_data/period,lw=5) #,'o--') #,markevery=100)
    #plt.xlim(155,157)
    ax2.set_xlabel("time")
    ax2.set_ylabel(r"y/$\lambda$")
    ax1.set_ylabel(r"F(t)")

    ax1.set_xticks([])

    ax1.set_xlim(0,time_data[-1]+1)
    ax2.set_xlim(0,time_data[-1]+1)
    #ax2.set_ylim(0,y_data/period[-1]+0.1)

    #add horizontal lines for potential minima
    for i in [1,3,5]:
        #ax2.hline(i)
        #, label = "U(y) = -U$_0$")
        ax2.axhline(y=i, color = "black", linestyle = "--")

    for i in [1,5,9,13]:
        #peaks of the FD curve
        ax2.axvline(x=i/(4*freq), color = "black", linestyle = "--") 
        ax1.axvline(x=i/(4*freq), color = "black", linestyle = "--") 

    ax1.text(0.02,0.9,"(a)",transform = ax1.transAxes,backgroundcolor="white")
    ax2.text(0.02,0.9,"(b)",transform = ax2.transAxes,backgroundcolor="white")
        
    plt.tight_layout(pad=0.2)
    plt.savefig(parameters['filename'])
    
    return 


#####################################################################
#def plot_velocity_force(time_data,FDC_data,vy_data,parameters):
#    '''
#    '''

#####################################################################
def plot_velocity_force(avg_FDC_data,avg_vy_data,parameters): #velocity_data,FDC_data,p):
    '''Make Fig. 3 in AJP
    '''

    #for plot
    #calculate the driving force from F_DC and time
    #F_AC = p['F_AC']
    #freq = p['freq']
    #period = p['period']
    
    fig = plt.figure(figsize=(6,4))
    gs=gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles

    #calculate from np.array rather than incrementally from subroutines
    #F_drive = FDC_data + F_AC*np.sin(2*np.pi*freq*time_data)

    #max occurs at
    # 2*np.pi*freq*time_data = m*pi/2, where m = 1, 5, 9, 13, ...
    # time = m/(4*freq)

    #plot the data
    #ax1.plot(FDC_data,vy_data) #,label="F$^{dc}$")
    ax1.plot(avg_FDC_data,avg_vy_data,'.') #,label="F$^{dc}$")

    #ax1.legend(loc=4,fontsize=20,borderaxespad=0.0,frameon=0,handlelength=1.5) #,labelspacing=0.2

    #plt.xlim(155,157)
    ax1.set_xlabel("F$^{dc}$")
    ax1.set_ylabel(r"v$_y$")

    #ax1.set_xlim(0,time_data[-1]+1)
    #ax2.set_ylim(0,y_data/period[-1]+0.1)
        
    plt.tight_layout(pad=0.2)
    plt.savefig(parameters['filename'])
    
    return 

#-------------------------------------------------------------------
def single_particle(parameters,plot="y-position"):
    '''Run MD simulation for a single particle driven across a washboard potential with an oscillating driving force

    Required argument(s): 
    p type: dictionary
    '''

    #define the constants (parameters and initial conditions)
    dt = parameters['dt']
    vy = parameters['vy0']
    avg_vy = 0 #vy
    time = parameters['time0']

    #start the particle in a local minima
    y = parameters['y0']

    F_DC = parameters['F_DC']
    F_DC_incr = parameters['F_DC_incr']
    F_DC_max = parameters['F_DC_max']
    
    #integer time steps
    maxtime=parameters['maxtime']
    writemovietime=parameters['writemovietime']
    
    #define arrays to hold data as a function of time
    array_length = int(maxtime/writemovietime)
    y_data = np.zeros(array_length)
    FDC_data = np.zeros(array_length)
    vy_data = np.zeros(array_length)
    time_data = np.zeros(array_length)

    avg_vy_data = np.zeros(int(maxtime/parameters['decifactor']))
    avg_FDC_data = np.zeros(int(maxtime/parameters['decifactor']))

    #loop through the integer time steps in the simulation
    for int_time in range(0,maxtime):

        #increase F_DC if needed
        if int_time > 0 and abs(F_DC - F_DC_max) > F_DC_incr: 
             F_DC = ramp_dc_force(F_DC,int_time, parameters)

        #apply the force calculations for the current position/time
        #note vy is not necessary 
        y, vy, avg_vy = md_step(y, time, F_DC, avg_vy, parameters)

        if int_time % writemovietime == 0:
            #print(int_time)
            #update the data
            i = int(int_time/writemovietime)
            y_data[i] = y
            vy_data[i] = vy
            FDC_data[i] = F_DC
            time_data[i] = time

        if int_time % parameters['decifactor'] == 0:
            j = int(int_time/parameters['decifactor'])

            #print(int_time,len(avg_vy_data))
            avg_vy_data[j] = avg_vy
            avg_FDC_data[j] = F_DC

            #reset the running average
            avg_vy = 0

        #update the time (simulation units)
        time += dt

    if plot == "y-position":
        plot_force_position_vs_time(time_data,FDC_data,y_data,parameters)
        return
    
    elif plot == "y-velocity":  
        return np.average(avg_vy_data)
      



    
#########################################################################
def set_parameters():
    '''Hardcode simulation parameters here.  In a compiled language like C this would be in a separate file that would be read into the code at execution time, so you wouldn't need to recompile to change the parameters
    '''

    #use a python dict (i.e. a hash to organize the parameters)
    dict={}

    dict['dt'] = 0.001 #timestep in simulation units

    dict['Sy']=36.5         #height of system
    dict['Sx']=36.5         #width of system
    dict['vy0'] = 0.0       #initial velocity of particle
    dict['time0'] = 0       #initial time in integer timesteps

    #start the particle in a local minima (not a local max!)
    dict['y0'] = 1. #period  #Sy/2

    #landscape potential
    dict['F0'] = 0.1
    dict['Np'] = 20         #number of troughs in the substrate
    dict['period'] = dict['Sy']/dict['Np']  #spatial period of substrate in y-direction

    #control the "constant" component of driving force
    dict['F_DC'] = 0              #initial F_DC
    dict['F_DC_incr'] = 0.01      #amount to increase FDC at every drop step
    dict['F_DC_max'] = 0.1        #"constant" driving force for most of simulation
    dict['drop'] = 4000           #integer timesteps to "ramp" the DC force
    dict['decifactor'] = 4000     #decimation factor integer timesteps to average force

    #control the oscillating component of driving force
    dict['F_AC'] = 0.05              #amplitude of force oscillation
    dict['freq'] = 0.01              #frequence of force osillation
    
    #integer time steps
    dict['maxtime']=400000        #total time steps in simulation
    dict['writemovietime']=1000   #interval to write data to arrays for plotting
    return dict


    
#-------------------------------------------------------------------
if __name__ == "__main__":

    parameters = set_parameters()

    make_fig2 = False
    make_fig3 = True

    #--------------------------------------------------------------
    #make Fig. 2 - run a single particle at a single driving force
    #--------------------------------------------------------------

    if make_fig2:
        parameters['filename']="single_particle.pdf"

        #run a single MD simulation for a set of parameters
        single_particle(parameters)

    if make_fig3: 
        #reset some parameters for subsequent figures
        #parameters['F_DC_incr'] = 0.01        #amount to increase FDC at every drop step
        #parameters['drop'] = 2000            #integer timesteps to "ramp" the DC force
        #parameters['decifactor'] = 2000      #integer timesteps to "ramp" the DC force

        #parameters['maxtime']=3000000        #total time steps in simulation
        parameters['filename']="sweep_FDC_vs_vx.pdf"
        
        delta_Fdc = 0.001
        Fdc_max=0.4
        max_value = int(Fdc_max/delta_Fdc)
        avg_vy_data = np.zeros(max_value)
        Fdc_data = np.arange(0,Fdc_max,delta_Fdc)

        #print(Fdc_data)
        #sys.exit()

        for i in range(len(Fdc_data)):

            #reset - though I think unnecessary
            parameters['F_DC'] = 0              #initial F_DC

            #"constant" driving force for most of simulation
            parameters['F_DC_max'] = Fdc_data[i]

            #parameters['y0'] = 5.0*parameters['period']

            #run a single MD simulation for a set of parameters
            avg_vy_data[i] = single_particle(parameters,plot="y-velocity")
            

        plot_velocity_force(Fdc_data,avg_vy_data,parameters)
    
    sys.exit()

