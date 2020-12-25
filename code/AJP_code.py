#################################################################
'''
Molecular dynamics simulation of synchronization in driven particles

Danielle McDermott
Tiare Guerrero
'''
#################################################################
import numpy as np
import math
import matplotlib.pyplot as plt

##################################################################
def landscape_potential(V0,y,period):
    '''
    '''
    
    return V0*np.cos(2*np.pi*y/period)

###################################################################
def external_drive(F_AC,F_DC,frequency,time):
    '''
    '''
    
    return F_DC + F_AC*np.sin(2*np.pi*frequency*time)

#####################################################################
def md_step(y, vy, dt, time):
    '''
    '''

    F_AC = 0.05
    F_DC = 0.05
    V0 = 0.5
    period = 1.0
    frequency = 0.1
    
    #integrate the equation of motion \eta v = F^net (eta=1)

    #substrate
    vy += landscape_potential(V0,y,period)

    #driving force
    vy += external_drive(F_AC,F_DC,frequency,time)

    #calculate the new position
    y += vy*dt
    
    return y, vy


#-------------------------------------------------------------------
if __name__ == "__main__":

    #define the constants (parameters and initial conditions)
    dt = 0.01 #hmm... not as in manuscript... single particle in a smoothly varying environment
    y = 0
    vy = 0
    time = 0.0

    #integer time steps
    maxtime=100

    #define arrays to hold data as a function of time
    y_data = np.zeros(maxtime)
    vy_data = np.zeros(maxtime)
    time_data = np.zeros(maxtime)

    #loop through the integer time steps in the simulation
    for int_time in range(0,maxtime):

        #apply the force calculations for the current position/time
        y, vy = md_step(y, vy, dt, time)

        #update the data
        y_data[int_time] = y
        vy_data[int_time] = vy
        time_data[int_time] = time

        #update the time (simulation units)
        time += dt


    plt.plot(time_data,y_data)
    plt.savefig("AJP.png")
