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
    
    #integrate the equation of motion \eta v = F^net

    vy += landscape_potential(V0,y,period)
    vy += external_drive(F_AC,F_DC,frequency,time)

    y += vy*dt
    
    return y, vy



if __name__ == "__main__":

    dt = 0.01
    y = 0
    vy = 0
    time = 0.0

    maxtime=100

    y_data = np.zeros(maxtime)
    vy_data = np.zeros(maxtime)
    time_data = np.zeros(maxtime)

    for int_time in range(0,maxtime):

        y, vy = md_step(y, vy, dt, time)

        y_data[int_time] = y
        vy_data[int_time] = vy
        time_data[int_time] = time

        time += dt


    plt.plot(time_data,y_data)
    plt.savefig("AJP.png")
