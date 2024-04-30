# Current state
# d/dt(N) + d/ds (g(s)N) = mu N
# mu(s) = 0
# g(s) = exp(-s/10)


import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def UPW_SPM(step, time, ds, dt, order, c):

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


    if order == 1:
        for t in range(0, len(time) - 1):
            # Time Splitting
            Ntemp = np.zeros([len(step)])

            # step 1 -- half time step
            for a in range(1,len(step)): 
                # Ntemp[a] = N[t,a] - (dt/(2*ds)) * (N[t,a] - N[t,a-1])
                Ntemp[a] = N[t,a] - (dt/(2*ds)) * (N[t,a] - N[t,a-1])

            # step 2 -- time step 
            for a in range(1,len(step)): 
                Ntemp[a] = Ntemp[a] * np.exp( - dt * mu[a] )      # exact solution 
                # Ntemp[a] = Ntemp[a] - mu[a] * Ntemp[a]            # numerical

            # step 3 -- half time step
            for a in range(1,len(step)): 
                N[t+1,a] = Ntemp[a] - (dt/(2*ds)) * (Ntemp[a] - Ntemp[a-1])
            
            # without time-splitting        
            # for a in range(1,len(step)):         
                # N[t+1,a] = N[t,a] - dt *( ( (N[t,a] - N[t,a-1]) / ds ) + mu[a] * N[t,a] ) # first order upwind method
   
    dtda = dt/ds 

    if order == 2:

        for t in range(0, len(time) - 1):

            exact_sol = np.exp( - dt * mu ) 
            
            # Time Splitting
            Ntemp = np.zeros([len(step)])

            # step 1 -- half time step
            for a in range(2,len(step)): 
                dNda = (3 * N[t,a] - 4 * N[t,a-1] + N[t,a-2]) 
                Ntemp[a] = N[t,a] - (dt/(4*ds)) * dNda 

            # step 2 -- time step 
            Ntemp[:] = Ntemp[:] * exact_sol[:]      # exact solution 

            # step 3 -- half time step
            for a in range(2,len(step)): 
                dNda = (3 * Ntemp[a] - 4 * Ntemp[a-1] + Ntemp[a-2]) 
                N[t+1,a] = Ntemp[a] - (dt/(4*ds)) * dNda 


    # trying to make it faster
    # if order == 2:
    #     # Precompute constants
    #     quarter_dt_da = dt / (4 * ds)
    #     exp_mu_dt = np.exp(-dt * mu)

    #     for t in range(len(time) - 1):
    #         # Time Splitting
    #         Ntemp = np.zeros_like(N[t])

    #         # Step 1 -- half time step
    #         for a in range(2, len(step)):
    #             dNda = 3 * N[t, a] - 4 * N[t, a - 1] + N[t, a - 2]
    #             Ntemp[a] = N[t, a] - quarter_dt_da * dNda

    #         # Step 2 -- time step
    #         Ntemp[2:] *= exp_mu_dt[2:]

    #         # Step 3 -- half time step
    #         for a in range(2, len(step)):
    #             dNda = 3 * Ntemp[a] - 4 * Ntemp[a - 1] + Ntemp[a - 2]
    #             N[t + 1, a] = Ntemp[a] - quarter_dt_da * dNda


         
        # # Second-order upwind scheme to solve the PDE
        # for t in range(len(time) - 1):
        #     for i in range(2, len(step)):
        #         # Second-order upwind difference for spatial derivative
        #         dNda = (3*N[t, i] - 4*N[t, i-1] + N[t, i-2]) / (2 * ds)
        #         # Time derivative and decay
        #         dNdt = -dNda - mu[i] * N[t, i]
        #         # Update using Euler's method
        #         N[t+1, i] = N[t, i] + dt * dNdt

        #     # Boundary conditions, assuming no flux at the boundary
        #     N[t+1, 0] = N[t+1, 1]  # Reflective or no-flux boundary condition at a_min
        #     N[t+1, 1] = N[t+1, 2]  # This can be adjusted based on physical assumptions


        # for t in range(0, len(time) - 2):
        #     if t == 0:
        #         print('using 1st order upwind')

        #         for a in range(1,len(step)):    
        #             N[t+1,a] = N[t,a] - dt *( ( (N[t,a] - N[t,a-1]) / ds ) + mu[a] * N[t,a] ) # first order upwind method
                
        #     else:    
        #         # print('using 2nd order upwind')
        #         for a in range(2,len(step)): 

        #             # first order time, second order step 
        #             dNda = (3 * N[t,a] - 4 * N[t,a-1] + N[t,a-2]) / (2 * ds)

        #             N[t+1,a] = N[t,a] - dt * (dNda + mu[a] * N[t,a])

    # if order == 2:
    #     for t in range(0, len(time) - 2):
    #         if t == 0:
    #             for a in range(1, len(step)):
    #                 N[t+1, a] = N[t, a] - dt * ((N[t, a] - N[t, a-1]) / ds + mu[a] * N[t, a])
    #         else:
    #             for a in range(2, len(step)):
    #                 # Apply the second-order upwind discretization
    #                 advection_term = (3 * N[t, a] - 4 * N[t, a-1] + N[t, a-2]) / (2 * ds)
    #                 N[t+1, a] = N[t, a] - dt * (advection_term + mu[a] * N[t, a])
                    
    return N