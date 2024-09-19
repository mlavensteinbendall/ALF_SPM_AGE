import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def LW_SPM(step, time, ds, dt, c):

    # inital condition -- population at t=0
    N = np.zeros([len(time),len(step)])
    N[0,:] = np.exp(-(step - 5)**2) 

    # mortality rate
    if c == 0:
        mu = np.zeros([len(step)])
    elif c == 1:
        mu = np.ones([len(step)])
    else:
        for i in range(len(step)):
            mu[i] = c


    for t in range(0, len(time)-1):

        # Time Splitting
        Ntemp  = np.zeros([len(step)])
        Ntemp2 = np.zeros([len(step)])

        # step 1 -- half time step
        for s in range(1,len(step)-1): 
            first_centeral_diff  = (N[t,s+1]              - N[t,s-1]) / (2*ds)
            second_centeral_diff = (N[t,s+1] - 2 * N[t,s] + N[t,s-1]) / ds**2

            Ntemp[s] = N[t,s] - (dt/2) * first_centeral_diff + (dt**2/8) * second_centeral_diff
            # Ntemp[a] = N[t,a] - (dt/(2*ds)) * (N[t,a] - N[t,a-1]) + (dt**2/(2**3 * ds**2)) * (N[t,a] - 2 * N[t,a-1] + N[t,a-2])  

        # step 2 -- full time step 
        for s in range(1,len(step)-1): 
            Ntemp2[s] = Ntemp[s] * np.exp( - dt * mu[s] )      # exact solution 

        # step 3 -- half time step
        for s in range(1,len(step)-1):
            # first_centeral_diff = (N[t,s+1] - N[t,s-1]) / (2*ds)                # update this to 
            # second_centeral_diff = (N[t,s+1] - 2 * N[t,s] + N[t,s-1]) / ds**2   # update this to 

            first_centeral_diff  = (Ntemp2[s+1]                 - Ntemp2[s-1]) / (2*ds)  # updated 
            second_centeral_diff = (Ntemp2[s+1] - 2 * Ntemp2[s] + Ntemp2[s-1]) / ds**2   # updated 

            N[t+1, s] = Ntemp2[s] - (dt/2) * first_centeral_diff + (dt**2/8) * second_centeral_diff

        # Set boundaries to zero
        N[t + 1, 0] = 0
        N[t + 1, -1] = 0
                
    return N




        # Ntemp[0] = np.exp(-(step[0]-5)**2)
        # Ntemp2[0] = np.exp(-(step[0]-5)**2)

        # does not work
        # step 1 -- half time step
        # for a in range(0,len(step)): 
        #     Ntemp[a] = N[t,a] * np.exp( -mu[a] * dt * 0.5 ) 
        #     # print(np.exp( -mu[a] * dt * 0.5 ) ) gives 1
        #     # print(Ntemp)
        #     # print(N[t,a])

        # # step 2 -- time step
        # for a in range(1, len(step) - 1): 
        #     first_centeral_diff     = (Ntemp[a+1]                - Ntemp[a-1]) * (  0.5 / ds      )
        #     second_centeral_diff    = (Ntemp[a+1] - 2 * Ntemp[a] + Ntemp[a-1]) * (  1   / ds**2   )

        #     # Ntemp[a] = Ntemp[a] - (dt) * first_centeral_diff + (0.5 * dt**2) * second_centeral_diff

        #     # first_centeral_diff     = (N[t,a+1]                - N[t,a-1]) * (  0.5 / ds      )
        #     # second_centeral_diff    = (N[t,a+1] - 2 * N[t,a] + N[t,a-1]) * (  1   / ds**2   )

        #     Ntemp[a] = Ntemp[a] - (dt) * first_centeral_diff + (0.5 * dt**2) * second_centeral_diff


        # # step 3 -- half time step
        # for a in range(0,len(step)): 
        #     N[t+1, a] = Ntemp[a] * np.exp( -mu[a] * dt * 0.5 ) 