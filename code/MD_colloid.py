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

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#set default font of plot
plt.rc('font', size=24)


##################################################################
def brownian_distribution(parameters):
    '''calculates a series of normally distributed random numbers 

    Required Parameters
    parameters: dict
    contains all constants of the simulation, notably the substrate period and barrier height AP
    '''
    mu, sigma = 0.0, 1.0 # mean and standard deviation

    #scale the prefactor by 2 \eta kb T in simulation units    
    #where 'temperature' has units of force
    kb_T = parameters['temperature']

    #one force for every timestep
    timesteps = int(parameters['maxtime'])
    
    #use numpy's built-in function to generate a randomized force
    #do so in bulk for every timestep in the simulation
    fy_rand = kb_T * np.random.normal(mu, sigma, timesteps)

    return fy_rand
    
    
##################################################################
def landscape_force(y,parameters):
    '''landscape potential is equation (1) in the manuscript, 
    U(y) = U0 cos(Np*2pi*y/L)
    while force is equation (5) - for now anyway
    F_y  = -dU/dy = -U0/period * sin(Np*2pi*y/L)
    potential minima occur at U(y) = -U0, or cos(Np*2pi*y/L) = -1
    meaning Np*2pi*y/L = m pi, where m is odd
    
    Required Parameters
    y (float): location of particle
    parameters (dict): all constants of the simulation, notably the substrate period and barrier height AP

    Optional arguments
    none
    '''

    period = parameters['period']
    AP = parameters['AP']
    
    return -AP*np.sin(2*np.pi*y/period)

###################################################################
def external_drive(time, parameters):
    '''Apply the driving force

    Required Parameters
    time (float): 
    parameters (dict): all constants of the simulation

    Optional arguments
    none
    '''

    F_AC = parameters['F_AC']
    F_DC = parameters['F_DC']
    frequency = parameters['freq']
    
    F_ext = F_DC + F_AC*np.sin(2*np.pi*frequency*time) #
    
    return F_ext


#####################################################################
def average_velocity(int_time, vy, avg_vy, parameters):
    '''calculate the average velocity over 
    a range of N timesteps, where N=decifactor 
    Use all vy values, not just those 
    saved during "writemovietime.  
    '''

    decifactor = parameters['decifactor']

    avg_vy += vy
            
    if int_time % parameters['decifactor'] == 0: 
        avg_vy /= decifactor
        #print(int_time, avg_vy)
        
    return avg_vy


#####################################################################
def md_step(y, int_time, avg_vy, parameters, ft=0):
    '''
    Required Arguments:
    y (float): position of particle
    int_time (int): count of MD steps
    avg_vy (float): calculated value (why/how?)
    parameters (dict): all the things

    Optional Arguments:
    ft (float) random kick due to temperature
    '''

    dt = parameters['dt']
    SY = parameters['Sy']
    F_DC = parameters['F_DC']
    
    temperature = parameters['temperature']

    time = int_time * dt
    
    #integrate the equation of motion \eta v = F^net (eta=1)

    #reset vy for every timestep since ay = 0 
    vy = landscape_force(y, parameters)

    #driving force
    vy += external_drive(time, parameters)

    if ft != 0:
        vy += ft

    #calculate the average velocity over all vy values
    #if decifactor > 1, we are accumulating many values to
    #create a smooth curve and smaller dataset
    if parameters['decifactor'] > 1:
        avg_vy = average_velocity(int_time, vy, avg_vy, parameters)
    else:
        avg_vy = vy

    #calculate the new position
    y += vy*dt

    #check periodic boundary conditions
    #if we've left the edge of the system, remap 
    if y > SY:
        y -= SY
    elif y < 0:
        y += SY
    
    return y, vy, avg_vy


#####################################################################
def plot_position_vs_time(ax2,time_data,y_data,p):
    '''Brownian figure, Fig. 5
    '''

    #for plot
    #calculate the driving force from F_DC and time
    F_AC = p['F_AC']
    freq = p['freq']
    period = p['period']

    ax2.plot(time_data,y_data/period,lw=5) #,'o--') #,markevery=100)
    #plt.xlim(155,157)
    ax2.set_xlabel(r"time ($\tau$)")
    ax2.set_ylabel(r"y/$\lambda$")

    ax2.set_xlim(0,time_data[-1]+1)
    #ax2.set_ylim(-0.1,y_data[-1]/period+0.2)

    #ax2.text(0.9,0.84,"(b)",transform = ax2.transAxes,backgroundcolor="white",zorder=-10)
    
    #add horizontal lines for potential minima

    '''
    for i in [1,3,5,7]:
        ax2.axhline(y=i, color = "black", linestyle = "--")

    for i in [1,5,9,13]:
        ax2.axvline(x=i/(4*freq), color = "black", linestyle = "--") 

    ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(8))
    ax1.grid(axis='both',which='both')
    '''
    #ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
    #ax2.xaxis.set_minor_locator(AutoMinorLocator(8))
    #ax2.grid(axis='both',which='both')
    
    return 


#####################################################################
def plot_force_position_vs_time(time_data,FDC_data,y_data,p):
    '''Make Fig. 2 

    panel (a) is driving force vs. time
    panel (b) is y/period vs. time
    '''

    #for plot
    #calculate the driving force from F_DC and time
    F_AC = p['F_AC']
    freq = p['freq']
    period = p['period']
    
    fig = plt.figure(figsize=(8,8)) #.25))
    gs=gridspec.GridSpec(2,1)
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
    ax2 = fig.add_subplot(gs[1,0]) #,sharex=ax1)  #scatter plot of particles

    #calculate from np.array rather than incrementally from subroutines
    F_drive = FDC_data + F_AC*np.sin(2*np.pi*freq*time_data)

    #max occurs at
    # 2*np.pi*freq*time_data = m*pi/2, where m = 1, 5, 9, 13, ...
    # time = m/(4*freq)

    #plot the force vs. time data
    #ax1.plot(time_data,FDC_data,label="F$^{dc}$")
    ax1.plot(time_data,F_drive,label="F$^d$(t)",lw=5)
    
    #plot normalized position by period vs. time
    ax2.plot(time_data,y_data/period,lw=5) #,'o--') #,markevery=100)

    #set all the labels, ax1 and ax2 have the same time axes
    ax1.set_ylabel(r"$F^d(t)$")
    ax2.set_xlabel(r"time ($\tau$)")
    ax2.set_ylabel(r"y/$\lambda$")

    #ax1.set_xticklabels([])
    #ax1.xaxis.set_ticklabels([])
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    ax1.set_xlim(0,time_data[-1]+1)
    ax2.set_xlim(0,time_data[-1]+1)
    #ax2.set_ylim(-0.1,y_data[-1]/period+0.2)

    ax1.text(0.9,0.8,"(a)",transform = ax1.transAxes,backgroundcolor="white")
    ax2.text(0.9,0.06,"(b)",transform = ax2.transAxes,backgroundcolor="white")

    #ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    #ax2.xaxis.set_minor_locator(AutoMinorLocator(1))
    #ax1.xaxis.set_minor_locator(AutoMinorLocator(1))
    ax2.grid(axis='both',which='both')
    ax1.grid(axis='both',which='both')

        
    plt.tight_layout(pad=0.1) #,h_pad=-0.3)
    plt.savefig(parameters['filename'])
    
    return 


#####################################################################

def plot_phase(ax,y_data,avg_vy_data,p):
    '''Make Fig. 8 in AJP
    '''

    #for plot
    #calculate the driving force from F_DC and time

    #plot the position within the period
    #y_data /= p['period']
    vbar = np.average(avg_vy_data)
    time = np.arange(1,p['maxtime']+1,1)*p['dt']
        
    plot_delay = 100   
    #no subtraction
    
    # driven phase variables
    phi_y_data = 2*np.pi*(y_data - vbar*time)/(p['period'])
    dot_phi_y_data = 2*np.pi*(avg_vy_data - vbar)/(p['period'])

    if 1:
        #plot the driven phase variables
        ax.scatter(phi_y_data[plot_delay:],dot_phi_y_data[plot_delay:],s=5) #,lw=5)
    else:
        #non-driven phase variables
        ax.scatter(y_data[plot_delay:],avg_vy_data[plot_delay:]) #,lw=5)

    ax.set_xlabel(r"$\phi_y(t)$")
    ax.set_ylabel(r"$\dot{\phi}_y(t)$")
    #ax.yaxis.set_label_coords(-0.2, 0.5)

    #ax1.set_xticks([])
    #ax1.set_xlim(0,time_data[-1]+1)
    
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
    #avg_FDC_data = np.zeros(int(maxtime/parameters['decifactor']))

    #loop through the integer time steps in the simulation
    for int_time in range(0,maxtime):

        #apply the force calculations for the current position/time
        #note vy is not necessary 
        y, vy, avg_vy = md_step(y, int_time, avg_vy, parameters)

        #collect up all the information to plot
        if int_time % writemovietime == 0:
            #print(int_time)
            #update the data
            i = int(int_time/writemovietime)
            y_data[i] = y
            vy_data[i] = vy
            FDC_data[i] = F_DC
            time_data[i] = time

        #if we're at the decimation factor
        if int_time % parameters['decifactor'] == 0:

            #find the integer for the array
            j = int(int_time/parameters['decifactor'])

            #save the running average 
            avg_vy_data[j] = avg_vy

            #this isn't the average, merely the value.  do i 
            #avg_FDC_data[j] = F_DC

            #reset the running average
            avg_vy = 0

        #update the time (simulation units)
        time += dt

    if plot == "y-position":
        plot_force_position_vs_time(time_data,FDC_data,y_data,parameters)
        return

    elif plot == "phase":
        #plot_phase(y_data,avg_vy_data,parameters)
        return y_data, avg_vy_data
    
    elif plot == "just_position":
        plot_position_vs_time(parameters['axis'], time_data,y_data,parameters)
        return
    
    elif plot == "y-velocity":
        #print(np.average(avg_vy_data), avg_vy_data)
        return np.average(avg_vy_data) #/decifactor
      

#----------------------------------------------------------------
#
#----------------------------------------------------------------
def brownian_particle(parameters,ax):
    '''
    '''

    #define the constants (parameters and initial conditions)
    dt = parameters['dt']
    vy = parameters['vy0']
    avg_vy = 0 #vy
    time = parameters['time0']

    #start the particle in a local minima
    y = parameters['y0']

    #assume no applied force for now
    F_DC = 0.0 #parameters['F_DC']
    parameters['F_AC']=0.0
    
    #integer time steps
    maxtime=parameters['maxtime']
    writemovietime=parameters['writemovietime']
    
    #define arrays to hold data as a function of time
    array_length = int(maxtime/writemovietime)
    y_data = np.zeros(array_length)
    #FDC_data = np.zeros(array_length)
    vy_data = np.zeros(array_length)
    time_data = np.zeros(array_length)

    avg_vy_data = np.zeros(int(maxtime/parameters['decifactor']))
    #avg_FDC_data = np.zeros(int(maxtime/parameters['decifactor']))

    #array should be the same length as above
    fy_rand = brownian_distribution(parameters)
    
    #loop through the integer time steps in the simulation
    for int_time in range(0,maxtime):


        #apply the force calculations for the current position/time
        #note vy is not necessary

        y, vy, avg_vy = md_step(y, int_time,
                                avg_vy, parameters, ft=fy_rand[int_time])
        
        if int_time % writemovietime == 0:
            #print(int_time)
            #update the data
            i = int(int_time/writemovietime)
            y_data[i] = y
            vy_data[i] = vy
            #FDC_data[i] = F_DC
            time_data[i] = time

        if int_time % parameters['decifactor'] == 0:
            j = int(int_time/parameters['decifactor'])

            #print(int_time,len(avg_vy_data))
            avg_vy_data[j] = avg_vy

            #print(int_time,avg_vy_data[j])

            #avg_FDC_data[j] = F_DC

            #reset the running average
            avg_vy = 0

        #update the time (simulation units)
        time += dt

    plot_position_vs_time(ax,time_data,y_data,parameters)
    return
    
    
#########################################################################
def set_parameters():
    '''Hardcode simulation parameters here.  In a compiled language like C this would be in a separate file that would be read into the code at execution time, so you wouldn't need to recompile to change the parameters
    '''

    #use a python dict (i.e. a hash to organize the parameters)
    dict={}

    dict['dt'] = 0.1 #005 #timestep in simulation units
    
    #integer time steps
    dict['maxtime']=4000        #total time steps in simulation
    dict['writemovietime']=1   #interval to write data to arrays for plotting

    dict['Sy']=36.5         #height of system
    dict['Sx']=36.5         #width of system
    dict['vy0'] = 0.0       #initial velocity of particle
    dict['time0'] = 0       #initial time in integer timesteps

    #start the particle in a local minima (not a local max!)
    dict['y0'] =  0 #dict['Sy']/2

    #landscape potential - Fig2 - 0.1
    dict['AP'] = 0.1 #
    dict['Np'] = 20         #number of troughs in the substrate
    dict['period'] = dict['Sy']/dict['Np']  #spatial period of substrate in y-direction

    #Brownian motion factor
    dict['temperature'] = 5.7*dict['AP']
    
    dict['F_DC'] = 0.05 #0.1
    dict['drop'] = dict['maxtime']   #integer timesteps to "ramp" the DC force
        
    dict['decifactor'] = 1 #5000     #decimation factor integer timesteps to average force

    #control the oscillating component of driving force

    #fig2
    dict['F_AC'] = 0.07 #0.1              #amplitude of force oscillation
    dict['freq'] = 0.01              #frequence of force osillation
    
    #parameters
    #dict['F_AC'] = 0.8              #amplitude of force oscillation
    #dict['freq'] = 0.006              #frequence of force osillation
    

    return dict

#-------------------------------------------------------------------
if __name__ == "__main__":

    parameters = set_parameters()

    #select which figure in the AJP you would like to make
    #figures 2 to 8
    make_fig = 3

    parameters['filename']="fig%d.pdf"%(make_fig)

    letter=["(a)","(b)","(c)","(d)","(e)","(f)","(g)"]

    #--------------------------------------------------------------
    #run a single particle at a single driving force
    #--------------------------------------------------------------
    if make_fig == 2 or make_fig == 8:

        #run a single MD simulation for a set of parameters
        if make_fig == 2:
            single_particle(parameters)
            
        elif make_fig == 8:
            #phase figure!!!!
            #parameters['drop'] = 4000   
            parameters['maxtime'] = 3000


            #parameters['F_AC'] = 0.2
            #parameters['freq'] = 0.01
            #parameters['AP'] = 0.05
            parameters['y0'] = 0 #parameters['period']/2

            fig = plt.figure(figsize=(10,12))
            gs=gridspec.GridSpec(3,2)

            FDC = [0.05, 0.1,0.12, 0.15, 0.2, 0.3]
            k=0
            for i in range(3):
                for j in range(2):
                    ax = fig.add_subplot(gs[i,j])  #scatter plot of particles
                    ax.text(0.02,0.9,letter[k],
                            transform = ax.transAxes,
                            backgroundcolor="white",zorder=-10)

                    parameters['F_DC'] = FDC[k]
                    #if k == 5:
                    #    parameters['y0'] = parameters['period']
                
                    y_data, avg_vy_data = single_particle(parameters,plot="phase")
                    plot_phase(ax, y_data,avg_vy_data,parameters)
                    k+=1

            
            plt.tight_layout(pad=0.1,h_pad=-0.3) #,v_pad=-0.3)
            plt.savefig(parameters['filename'])

    #--------------------------------------------------------------
    #make shapiro steps with a range of F^{dc}
    #this is the place to include the strictly DC force...
    #--------------------------------------------------------------
    if make_fig == 3 or make_fig == 7: 

        #parameters['maxtime']=10000

        #RULZ = watz bads
        #if maxtime/decifactor not integer
        #if decifactor too bigz
        #parameters['decifactor'] = 2500 #parameters['maxtime']

        delta_Fdc = 0.001
        
        if make_fig == 3:

            Fdc_max=0.6+delta_Fdc
            parameters['freq']=0.1

            F_AC = [0.0, 0.07, 0.2]

        
        elif make_fig == 7:

            #IN EXPLORING THE HIGH FREQUENCY REGIME,
            #I'M SEEING STEP HEIGHTS LIKE JUNIPER (but still fractional)  
            parameters['freq'] = 0.1
            delta_Fdc = 0.01
            F_AC = [0.0,0.3]
            Fdc_max=0.5+delta_Fdc


        max_value = int(np.ceil(Fdc_max/delta_Fdc))
        avg_vy_data = np.zeros(max_value)
        Fdc_data = np.arange(0,Fdc_max,delta_Fdc)


        #make the figure 
        fig = plt.figure(figsize=(7,4))
        gs=gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles

        for F in F_AC:
            parameters['F_AC'] = F

            for i in range(len(avg_vy_data)):

                #sweep through a range of dc values
                parameters['F_DC'] = Fdc_data[i] 
                
                #run a single MD simulation for a set of parameters
                avg_vy_data[i] = single_particle(parameters,plot="y-velocity")            

            #normalize <vy> so you can count step number
            avg_vy_data /= (parameters['freq']*parameters['period'])
            print(avg_vy_data)
            
            #plot <vy> vs F_dc
            ax1.plot(Fdc_data,avg_vy_data,label=r"%1.2f"%(F))

            
            
        ax1.set_xlim(0,Fdc_data[-1])
        ax1.set_ylim(avg_vy_data[0]-0.001,avg_vy_data[-1])

        #plt.xlim(155,157)
        ax1.set_xlabel("F$^{dc}$")
        ax1.set_ylabel(r"$\langle$v$_y \rangle / \lambda f $")
        ax1.legend(loc='best',fontsize=20,numpoints=3,borderpad=0.1,handletextpad=0.2,title="F$^{ac}$")
                    
        plt.tight_layout(pad=0.2)
        plt.savefig(parameters['filename'])
    
    #--------------------------------------------------------------
    #study Brownian motion for increasing temperature
    #--------------------------------------------------------------
    if make_fig == 4:

        #place particle in center to reduce hops across PBC
        parameters['y0'] = parameters['Sy']/2

        #simulating for a long time, 
        parameters['dt'] = 0.1 #timestep in simulation units

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
            brownian_particle(parameters,ax)

            i+=1

        #save fig 4
        plt.tight_layout(pad=0.1,h_pad=-0.3)
        plt.savefig(parameters['filename'])
        
    #--------------------------------------------------------------
    #modify either frequency (fig5) or F^ac amplitude (fig6) or F^ac = 0
    #--------------------------------------------------------------
    if make_fig == 5 or make_fig == 6:
        
        parameters['dt'] = 0.1 #timestep in simulation units
        parameters['writemovietime']=1   #interval to write data to arrays for plotting


        if make_fig == 5:
            parameters['maxtime']=5000        #total time steps in simulation

            #we will sweep through the following independent variable

            frequency = [0.1, 0.05, 0.015, 0.005, 0.001]
            ind_var = frequency
            str_iv=" f=%1.3f"

        else:
            parameters['maxtime']=3000        #total time steps in simulation

            #we will sweep through the following independent variable
            F_AC = [0.2, 0.3, 0.4]
            ind_var = F_AC
            str_iv=" $F^{ac}$=%1.2f"
            parameters['freq'] = 0.01
           


        #make the figure
        n_plots = len(ind_var)
        fig = plt.figure(figsize=(10,3*n_plots))
        gs=gridspec.GridSpec(n_plots,1)
        
        #run a single MD simulation for a set of parameters in independent variable
        for i in range(n_plots):

            #make the subplot
            ax = fig.add_subplot(gs[i,0])
            ax.text(0.02,0.8,letter[i]+str_iv%(ind_var[i]),
                    transform = ax.transAxes,backgroundcolor="white")
            parameters['axis'] = ax

            #select the independent variable
            if make_fig == 5:
                parameters['freq'] = frequency[i]
            else:
                parameters['F_AC'] = F_AC[i]

                #horizontal grid so we can count hops
                ax.set_ylim(-0.9,20.9)
                ax.yaxis.set_minor_locator(AutoMinorLocator(5))
                ax.grid(axis='y',which='both')

            #run the MD simulation    
            single_particle(parameters,plot="just_position")

            if i < n_plots-1:
                ax.tick_params(labelbottom=False)
                ax.set_xlabel("")

        #save the figure
        plt.tight_layout(pad=0.1,h_pad=-0.3)
        plt.savefig(parameters['filename'])


    #--------------------------------------------------------------------------
    #animate this

    if make_fig == 10:

        print('animate')

    sys.exit()


####################################################################
def ramp_dc_force(F_DC, time, parameters):
    '''Not used
    '''
    F_DC_max = parameters['F_DC_max']
    drop = parameters['drop']
    F_DC_incr = parameters['F_DC_incr']
    
    if time % drop == 0 and F_DC < F_DC_max:
        F_DC += F_DC_incr

    return F_DC

#####################################################################
def plot_velocity_force(avg_FDC_data,avg_vy_data,parameters): #velocity_data,FDC_data,p):
    '''Make Fig. 3 in AJP
    '''
    
    ax1 = parameters['ax']
    ax1.plot(avg_FDC_data,avg_vy_data,'.') #,label="F$^{dc}$")


    return 

#sometimes convenient to start farther from the PBC
#parameters['y0'] = 5.0*parameters['period']


