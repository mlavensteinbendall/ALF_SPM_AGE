import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt


def LW_SPM(step, time, ds, dt, c):
    """Calculates the numerical solution using strang splitting, lax-wendroff, and runge-kutta method. 
    
    Args:
        step    (array): 
        time    (array):
        ds      (int):
        dt      (int):
        c       (int):
        
    Returns:
        N       (array of arrays): Represents the time
    """

    # inital condition -- population at t=0
    N = np.zeros([len(time),len(step)])
    N[0,:] = np.exp(-(step - 5)**2) 
    # print(N[0,:])

    # mortality rate
    if c == 0:
        mu = np.zeros([len(step)])
    elif c == 1:
        mu = np.ones([len(step)])
    else:
        for i in range(len(step)):
            mu[i] = c



## NUMERICAL SOLUTION 
    for t in range(0, len(time)-1):

        # Time Splitting
        Ntemp  = np.zeros([len(step)])
        Ntemp2 = np.zeros([len(step)])

        # step 1 -- half time step (Aging)
        for s in range(1,len(step)-1): 

            # Centeral Finite Difference
            first_centeral_diff  = (N[t,s+1]              - N[t,s-1]) / (2*ds)
            second_centeral_diff = (N[t,s+1] - 2 * N[t,s] + N[t,s-1]) / ds**2

            Ntemp[s] = N[t,s] - (dt/2) * first_centeral_diff + (dt**2/8) * second_centeral_diff

        # step 2 -- full time step (death)
        for s in range(0,len(step)): 
            # Ntemp2[s] = Ntemp[s] * np.exp( - dt * mu[s] )      # exact solution 

            # RK2
            k1        = Ntemp[s] * mu[s]                    # slope at current time
            N_star    = Ntemp[s] - dt * k1                  # intermediate value
            k2        = N_star * mu[s]                      # slope at intermediate value
            Ntemp2[s] = Ntemp[s] - (dt * 0.5) * (k1 + k2)   # update value

        # step 3 -- half time step (Aging)
        for s in range(1,len(step)-1):

            # Centeral Finite Difference
            first_centeral_diff  = (Ntemp2[s+1]                 - Ntemp2[s-1]) / (2*ds) # first order finite difference
            second_centeral_diff = (Ntemp2[s+1] - 2 * Ntemp2[s] + Ntemp2[s-1]) / ds**2  # second order finite difference

            N[t+1, s] = Ntemp2[s] - (dt/2) * first_centeral_diff + (dt**2/8) * second_centeral_diff

    # print("last time-step:", N[-1,:])

#------------------------

        # Set boundaries to zero
        N[t + 1, 0] = 0
        N[t + 1, -1] = 0
                
    return N


    # # slightly improved -- switched ordering
    # # for t in range(0, len(time)-1):

    # #     # Time Splitting
    # #     Ntemp  = np.zeros([len(step)])
    # #     Ntemp2 = np.zeros([len(step)])

    # #     # print("Ntemp shape = ", Ntemp.shape)

    # #     # step 1 -- half time step
    # #     for s in range(0,len(step)-1):
    # #         Ntemp[s] = N[t,s] * np.exp( - dt * 0.5 * mu[s] )      # exact solution 

    # #     # step 2 -- full time step 
    # #     for s in range(1,len(step)-1): 
    # #         first_centeral_diff  = (Ntemp[s+1]                 - Ntemp[s-1]) / (2*ds)  # updated 
    # #         second_centeral_diff = (Ntemp[s+1] - 2 * Ntemp[s] + Ntemp[s-1]) / ds**2   # updated 

    # #         Ntemp2[s] = Ntemp[s] - dt * first_centeral_diff + (dt**2/2) * second_centeral_diff

    # #     # step 3 -- half time step
    # #     for s in range(0,len(step)):
    # #         N[t+1, s] = Ntemp2[s] * np.exp( - dt * 0.5 * mu[s] )      # exact solution 

    #     # Set boundaries to zero
    #     N[t + 1, 0] = 0
    #     N[t + 1, -1] = 0
                
    # return N




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