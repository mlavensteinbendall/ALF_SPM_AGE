import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def LW_SPM(step, time, ds, dt, c):

    # mortality rate
    if c == 0:
        mu = np.zeros([len(step)])
    elif c == 1:
        mu = np.ones([len(step)])
    else:
        for i in range(len(step)):
            mu[i] = c

    # inital condition -- population at t=0
    N = np.zeros([len(time),len(step)])
    N[0,:] = np.exp(-(step[:]-5)**2)


    for t in range(0, len(time) - 1):

        # Time Splitting
        Ntemp = np.zeros([len(step)])

        # does not work
        # # step 1 -- half time step
        # for a in range(len(step)): 
        #     Ntemp[a] = N[t,a] * np.exp( -mu[a] * dt * 0.5 ) 

        # # step 2 -- time step
        # for a in range(1,len(step)-1): 
        #     first_centeral_diff     = (Ntemp[a+1]                - Ntemp[a-1]) / (2*ds)
        #     second_centeral_diff    = (Ntemp[a+1] - 2 * Ntemp[a] + Ntemp[a-1]) / (ds**2)

        #     Ntemp[a] = Ntemp[a] - (dt) * first_centeral_diff + (dt**2/2) * second_centeral_diff

        # # step 3 -- half time step
        # for a in range(len(step)): 
        #     N[t+1, a] = Ntemp[a] * np.exp( -mu[a] * dt * 0.5 ) 


        # step 1 -- half time step
        for a in range(1,len(step)-1): 
            first_centeral_diff = (N[t,a+1] - N[t,a-1]) / (2*ds)
            second_centeral_diff = (N[t,a+1] - 2 * N[t,a] + N[t,a-1]) / ds**2

            Ntemp[a] = N[t,a] - (dt/2) * first_centeral_diff + (dt**2/8) * second_centeral_diff

            # Ntemp[a] = N[t,a] - (dt/(2*ds)) * (N[t,a] - N[t,a-1]) + (dt**2/(2**3 * ds**2)) * (N[t,a] - 2 * N[t,a-1] + N[t,a-2])

        # step 2 -- time step 
        for a in range(1,len(step)-1): 
            Ntemp[a] = Ntemp[a] * np.exp( - dt * mu[a] )      # exact solution 

        # step 3 -- half time step
        for a in range(1,len(step)-1):
            first_centeral_diff = (N[t,a+1] - N[t,a-1]) / (2*ds)
            second_centeral_diff = (N[t,a+1] - 2 * N[t,a] + N[t,a-1]) / ds**2

            N[t+1, a] = Ntemp[a] - (dt/2) * first_centeral_diff + (dt**2/8) * second_centeral_diff

        # for a in range(0,2):
        #     N[t+1,a] = N[t,a] 

        # without time splitting
        # for a in range(1,len(step)-1): 
        #     first_centeral_diff = (N[t,a+1] - N[t,a-1]) / (2*ds)
        #     second_centeral_diff = (N[t,a+1] - 2 * N[t,a] + N[t,a-1]) / ds**2

        #     N[t+1,a] = N[t,a] - (dt) * first_centeral_diff + (dt**2/2) * second_centeral_diff

              
    return N