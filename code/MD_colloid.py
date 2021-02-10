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
    '''
    '''
    mu, sigma = 0.0, 1.0 # mean and standard deviation

    #temperature has units of force too, because of course it does
    kb_T = parameters['temperature']

    timesteps = int(parameters['maxtime'])
    
    #use numpy's built-in function to generate a randomized force
    #do so in bulk for every timestep in the simulation
    #scale the prefactor by 2 \eta kb T in simulation units    
    fy_rand = kb_T * np.random.normal(mu, sigma, timesteps)

    return fy_rand
    
    #fx_rand = np.random.normal(mu, sigma, int(N_part*(N_part-1)))
    
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
def average_velocity(int_time, vy, avg_vy, parameters):
    '''calculate the average velocity over 
    a range of time values for N particles.  
    Use all vy values, not just those 
    saved during "writemovietime.  So far just "
    written for the single particle case, 
    where the code is not so interesting.  Subroutine 
    generally sums all vxi and vyi values across n particles and averages
    '''

    decifactor = parameters['decifactor']

    avg_vy += vy
            
    if int_time % parameters['decifactor'] == 0: 
        avg_vy /= decifactor
        #print(int_time, avg_vy)
        
    return avg_vy


#####################################################################
def md_step(y, int_time, F_DC, avg_vy, parameters, ft=0):
    '''

    ft random kicks due to temperature affects
    '''

    dt = parameters['dt']
    SY = parameters['Sy']
    temperature = parameters['temperature']

    time = int_time * dt
    
    #integrate the equation of motion \eta v = F^net (eta=1)

    #reset vy for every timestep since ay = 0 
    vy = landscape_force(y, parameters)

    #driving force
    vy += external_drive(F_DC, time, parameters)

    if ft != 0:
        vy += ft

    #calculate the average velocity over all vy values
    avg_vy = average_velocity(int_time, vy, avg_vy, parameters)

    '''
    decifactor = parameters['decifactor']
    if int_time == decifactor:
        print(int_time, avg_vy)
    '''
    #calculate the new position
    y += vy*dt

    #check periodic boundary conditions
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
    '''
    
    return 


#####################################################################
def plot_force_position_vs_time(time_data,FDC_data,y_data,p):
    '''Make Fig. 2 in AJP
    '''

    #for plot
    #calculate the driving force from F_DC and time
    F_AC = p['F_AC']
    freq = p['freq']
    period = p['period']
    
    fig = plt.figure(figsize=(8,7.25))
    gs=gridspec.GridSpec(2,1)
    ax1 = fig.add_subplot(gs[0,0])  #scatter plot of particles
    ax2 = fig.add_subplot(gs[1,0]) #,sharex=ax1)  #scatter plot of particles

    #calculate from np.array rather than incrementally from subroutines
    F_drive = FDC_data + F_AC*np.sin(2*np.pi*freq*time_data)

    #max occurs at
    # 2*np.pi*freq*time_data = m*pi/2, where m = 1, 5, 9, 13, ...
    # time = m/(4*freq)

        #plot the data
    #ax1.plot(time_data,FDC_data,label="F$^{dc}$")
    ax1.plot(time_data,F_drive,label="F$^d$(t)",lw=5)
    #ax1.legend(loc=4,fontsize=20,borderaxespad=0.0,frameon=0,handlelength=1.5) #,labelspacing=0.2
    ax2.plot(time_data,y_data/period,lw=5) #,'o--') #,markevery=100)
    #plt.xlim(155,157)
    ax2.set_xlabel(r"time ($\tau$)")
    ax2.set_ylabel(r"y/$\lambda$")
    ax1.set_ylabel(r"F(t)")

    ax1.set_xticks([])

    ax1.set_xlim(0,time_data[-1]+1)
    ax2.set_xlim(0,time_data[-1]+1)
    #ax2.set_ylim(-0.1,y_data[-1]/period+0.2)

    ax1.text(0.9,0.86,"(a)",transform = ax1.transAxes,backgroundcolor="white")
    ax2.text(0.9,0.06,"(b)",transform = ax2.transAxes,backgroundcolor="white")

    ax2.set_ylim(0.0,8.501)
    ax2.set_yticks([1,3,5,7])
    
    #add horizontal lines for potential minima
    #'''
    for i in range(1,10)[::2]:
        #ax2.hline(i)
        #, label = "U(y) = -U$_0$")
        ax2.axhline(y=i, color = "black", linestyle = "--")
    #'''
    for i in [1,5,9,13]:
        #peaks of the FD curve
        ax2.axvline(x=i/(4*freq), color = "black", linestyle = "--") 
        ax1.axvline(x=i/(4*freq), color = "black", linestyle = "--") 

        
    plt.tight_layout(pad=0.1,h_pad=-0.3)
    #g.tight_layout(fig, rect=[0, 0, 1, 1], h_pad=p_val,w_pad=p_val,pad=0.0)
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
    
    fig = plt.figure(figsize=(7,4))
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
    ax1.set_ylabel(r"$\langle$v$_y \rangle$")

    ax1.set_xlim(0,avg_FDC_data[-1])
    ax1.set_ylim(avg_vy_data[0]-0.001,avg_vy_data[-1])
        
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
        y, vy, avg_vy = md_step(y, int_time, F_DC, avg_vy, parameters)
        
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

            #print(int_time,avg_vy_data[j])

            avg_FDC_data[j] = F_DC

            #reset the running average
            avg_vy = 0

        #update the time (simulation units)
        time += dt

    if plot == "y-position":
        plot_force_position_vs_time(time_data,FDC_data,y_data,parameters)
        return
    
    elif plot == "just_position":
        plot_position_vs_time(parameters['axis'], time_data,y_data,parameters)
        return
    
    elif plot == "y-velocity":  
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
    F_DC_incr = 0.0 #parameters['F_DC_incr']
    F_DC_max = 0.0 #parameters['F_DC_max']
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

        #increase F_DC if needed
        #if int_time > 0 and abs(F_DC - F_DC_max) > F_DC_incr: 
        #     F_DC = ramp_dc_force(F_DC,int_time, parameters)

        #apply the force calculations for the current position/time
        #note vy is not necessary

        y, vy, avg_vy = md_step(y, int_time,
                                F_DC, avg_vy, parameters, ft=fy_rand[int_time])
        
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

    dict['dt'] = 0.001 #timestep in simulation units

    dict['Sy']=36.5         #height of system
    dict['Sx']=36.5         #width of system
    dict['vy0'] = 0.0       #initial velocity of particle
    dict['time0'] = 0       #initial time in integer timesteps

    #start the particle in a local minima (not a local max!)
    dict['y0'] =  0 #dict['Sy']/2

    #landscape potential - Fig2 - 0.1
    dict['F0'] = 0.1 #parameters 0.2
    dict['Np'] = 20         #number of troughs in the substrate
    dict['period'] = dict['Sy']/dict['Np']  #spatial period of substrate in y-direction

    #Brownian motion factor
    dict['temperature'] = 5.7*dict['F0']
    
    #integer time steps
    dict['maxtime']=400000        #total time steps in simulation
    dict['writemovietime']=1000   #interval to write data to arrays for plotting
    
    #control the "constant" component of driving force
    dict['F_DC_max'] = 0.1        #"constant" driving force for most of simulation
    if 0:
        dict['F_DC'] = 0              #initial F_DC
        dict['F_DC_incr'] = 0.01      #amount to increase FDC at every drop step
        dict['F_DC_max'] = 0.1        #"constant" driving force for most of simulation
        dict['drop'] = 4000           #integer timesteps to "ramp" the DC force
    else:
        dict['F_DC'] = dict['F_DC_max']  #initial F_DC
        dict['F_DC_incr'] = 0.0          #amount to increase FDC at every drop step
        dict['drop'] = dict['maxtime']   #integer timesteps to "ramp" the DC force
        
    dict['decifactor'] = 5000     #decimation factor integer timesteps to average force

    #control the oscillating component of driving force

    #fig2
    dict['F_AC'] = 0.05              #amplitude of force oscillation
    dict['freq'] = 0.01              #frequence of force osillation
    
    #parameters
    #dict['F_AC'] = 0.8              #amplitude of force oscillation
    #dict['freq'] = 0.006              #frequence of force osillation
    

    return dict

#-------------------------------------------------------------------
if __name__ == "__main__":

    parameters = set_parameters()

    #select which figure in the AJP you would like to make
    '''
    make_fig2 = False #True #False #True 
    make_fig3 = True #
    make_fig4 = False #True
    make_fig5 = False #True
    make_fig6 = False #True
    make_fig7 = False #True
    '''

    #Coded to make figures 2 to 8
    make_fig = 7

    letter=["(a)","(b)","(c)","(d)","(e)","(f)","(g)"]

    #--------------------------------------------------------------
    #run a single particle at a single driving force
    #--------------------------------------------------------------
    if make_fig == 2:
        
        parameters['dt'] = 0.1 #timestep in simulation units
        parameters['maxtime']=4000        #total time steps in simulation
        parameters['writemovietime']=1   #interval to write data to arrays for plotting
        parameters['decifactor']=1   #interval to write data to arrays for plotting

        parameters['filename']="single_particle_dt0.1.pdf"

        #run a single MD simulation for a set of parameters
        single_particle(parameters)

    #--------------------------------------------------------------
    #make shapiro steps with a range of F^{dc}
    #--------------------------------------------------------------
    if make_fig == 3 or make_fig == 7: 
        #reset some parameters for subsequent figures
        #parameters['F_DC_incr'] = 0.01        #amount to increase FDC at every drop step
        #parameters['drop'] = 2000            #integer timesteps to "ramp" the DC force
        #parameters['decifactor'] = 2000      #integer timesteps to "ramp" the DC force

        parameters['maxtime']=100000          #shorter and I don't get steps

        
        if make_fig == 3:
            delta_Fdc = 0.001
            Fdc_max=0.3+delta_Fdc
            parameters['filename']="fig3_sweep_FDC_vs_vx.pdf"

        elif make_fig == 7:
            delta_Fdc = 0.01
            parameters['drop'] = 10000    #integer timesteps to "ramp" the DC force
            parameters['decifactor'] = 10000      #integer timesteps to "ramp" the DC force

            parameters['F_AC'] = 0.4
            Fdc_max=0.6+delta_Fdc
            parameters['filename']="fig7_sweep_FDC_vs_vx.pdf"
            
        max_value = int(np.ceil(Fdc_max/delta_Fdc))
        avg_vy_data = np.zeros(max_value)
        Fdc_data = np.arange(0,Fdc_max,delta_Fdc)

        #print(Fdc_data)
        #sys.exit()

        for i in range(len(avg_vy_data)):

            #reset - though I think unnecessary
            parameters['F_DC'] = Fdc_data[i] #0              #initial F_DC

            #"constant" driving force for most of simulation
            parameters['F_DC_max'] = Fdc_data[i]

            #parameters['y0'] = 5.0*parameters['period']

            #run a single MD simulation for a set of parameters
            avg_vy_data[i] = single_particle(parameters,plot="y-velocity")            

        plot_velocity_force(Fdc_data,avg_vy_data,parameters)


    #--------------------------------------------------------------
    #study Brownian motion for increasing temperature
    #--------------------------------------------------------------
    if make_fig == 4:

        #place particle in center to remove hops across periodic boundary conditions
        parameters['y0'] = parameters['Sy']/2

        #simulating for a long time, 
        parameters['dt'] = 0.1 #timestep in simulation units
        parameters['maxtime']=3000000        #total time 
        parameters['writemovietime']=10000   #interval to write data 
        parameters['decifactor']=1   #interval to write data to arrays for plotting

        #save the figure in
        parameters['filename']="brownian_particle_dt0.1.pdf"

        #make the figure
        fig = plt.figure(figsize=(10,12))
        gs=gridspec.GridSpec(3,1)
        i=0

        #sweep through temperatures
        for temp in [4.0, 5.0, 6.0]:

            #create and label the subplot
            ax = fig.add_subplot(gs[i,0])
            ax.text(0.02,0.89,letter[i],
                    transform = ax.transAxes,backgroundcolor="white",zorder=-10)

            #Brownian motion factor as a ratio of the Ap value
            parameters['temperature'] = temp*parameters['F0']

            #run a single MD simulation for a set of parameters,
            #will be plotted on subplot
            brownian_particle(parameters,ax)

            i+=1

        #save fig 4
        plt.tight_layout(pad=0.1,h_pad=-0.3)
        plt.savefig(parameters['filename'])
        
    #--------------------------------------------------------------
    #modify either frequency (fig5) or F^ac amplitude (fig6)
    #--------------------------------------------------------------
    if make_fig == 5 or make_fig == 6:
        
        parameters['dt'] = 0.1 #timestep in simulation units
        parameters['maxtime']=3000        #total time steps in simulation
        parameters['writemovietime']=1   #interval to write data to arrays for plotting
        parameters['decifactor']=1   #interval to write data to arrays for plotting

        if make_fig == 5:
            parameters['filename']="parameters_fig5.pdf"
            frequency = [0.1, 0.05, 0.01, 0.005, 0.001]

            #we will sweep through the following independent variable
            ind_var = frequency
            str_iv=" f=%1.3f"

        else:
            F_AC = [0.2, 0.3, 0.4]

            ind_var = F_AC
            str_iv=" $F^{ac}$=%1.2f"
            
            #we will sweep through the following independent variable
            parameters['filename']="parameters_fig6.pdf"
            parameters['freq'] = 0.01

            #parameters['F0'] = 0.05              
            #parameters['F_DC_max'] = 0.05        #"constant" driving force for most of simulation


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


    #-----------------------------------------------------------------------------------
    sys.exit()

