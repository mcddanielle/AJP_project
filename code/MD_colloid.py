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

##################################################################
def landscape_potential(V0,y,period):
    '''but what is the landscape force?
    '''
    
    return V0*np.sin(2*np.pi*y/period)

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

    F_AC = 0.2
    #F_DC = 0.005 #5 #increment slowly
    period = 1.825
    V0 = 0.1 #*2*np.pi/period
    frequency = 0.01
    SX = 36.5
    SY = 36.5
    
    #integrate the equation of motion \eta v = F^net (eta=1)

    #reset vy for every timestep!  
    vy = landscape_potential(V0,y,period)

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
    y = 0
    vy = 0
    time = 0.0
    period = 1.825

    F_DC = 0
    F_DC_incr = 0.0 #01
    drop = 4000
    
    #integer time steps
    maxtime=200000
    writemovietime=1000
    
    #define arrays to hold data as a function of time
    array_length = int(maxtime/writemovietime)
    y_data = np.zeros(array_length)
    vy_data = np.zeros(array_length)
    time_data = np.zeros(array_length)

    #loop through the integer time steps in the simulation
    for int_time in range(0,maxtime):

        if int_time > 0: 
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
            time_data[i] = time

        #update the time (simulation units)
        time += dt


    #plot the data
    plt.plot(time_data,y_data/period,'o--') #,markevery=100)
    #plt.xlim(155,157)
    plt.xlabel("time")
    plt.ylabel(r"y/$\lambda$")
    plt.savefig("AJP.png")
