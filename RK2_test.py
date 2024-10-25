import numpy                as np
import matplotlib.pyplot    as plt 



def RK2_function(step, time, ds, dt, mu):

    # inital condition -- population at t=0
    N = np.zeros([len(time),len(step)])
    N[0,:] = np.exp(-(step - 5)**2) 


    ## NUMERICAL SOLUTION 
    for t in range(0, len(time)-1):

        for s in range(0,len(step)-1):
            
            Ntemp  = np.zeros([len(step)])
            # N[t+1, s] = N[t, s] * np.exp( - dt * mu[s] )      # exact solution 


            # # # RK2
            # k1        = N[t, s] * mu[s]                    # slope at beginning of interval
            # N_star    = N[t, s] - dt * k1                  # slope at the midpoint of interval
            # k2        = N_star * (mu[s] + mu[s+1])/2       # slope at end of interval
            # N[t+1, s] = N[t, s] - (dt * 0.5) * (k1 + k2)   # update value
            # # # N[t+1, s] = N[t, s] - (dt * k1) * (1 + dt/2)   # update value

            # # RK2
            # k1        = N[t, s] * mu[s]                     # slope at beginning of interval
            # N_star    = N[t, s] - dt * k1                  # slope at the midpoint of interval
            # k2        = N_star   * mu[s]                    # slope at end of interval
            # N[t+1, s] = N[t, s] - (dt * 0.5) * (k1 + k2)   # update value

            # RK2 for reaction term with mu(a) as a function of age CURRENT ONE
            k1        = N[t, s] * mu[s]                       # slope at the beginning of the time step (age = s)
            N_star    = N[t, s] - (dt / 2) * k1                # estimate N at the midpoint of the time step
            k2        = N_star * mu[s]                         # slope at the midpoint (age still = s)
            N[t+1,s]  = N[t,s] - dt * k2                      # update N for the full time step
            
            # N[t+1, s] = N[t, s] / (1 + dt * mu[s])
            # Crank-Nicolson update
            # N[t+1, s] = N[t, s] * (1 - (dt / 2) * mu[s]) / (1 + (dt / 2) * mu[s])

            # print(Ntemp)

            # N[t+1,s] = N[t, s] * np.exp( - mu[s] * time[t])
            # print(Ntemp)

            # N[t+1,s] = Ntemp[s]

            # LW
            # N[t+1, s] = N[t,s] - mu[s] * N[t,s] * dt + mu[s]**2 * N[t,s] * dt**2 / 2

            # Crank-Nicolson update
            # N[t+1, s] = N[t,s] * (1 - (dt / 2) * mu[s]) / (1 + (dt / 2) * mu[s])


            

        # Set boundaries to zero
        N[t + 1, 0] = 0
        N[t + 1, -1] = 0
                
    return N


# ds = 0.5
# dt = 0.01
# Smax = 15
# Tmax = 1

# time = np.arange(0, Tmax, dt)
# step = np.arange(0, Smax, ds)

# # inital condition -- population at t=0
# N = np.zeros([len(time),len(step)])
# N[0,:] = np.exp(-(step - 5)**2) 

# # mu = np.ones(len(step))
# mu = 0.5 * step


# for t in range(0, len(time) - 1):

#     # Time step (death)
#     for s in range(0,len(step)): 
#         # Ntemp2[s] = Ntemp[s] * np.exp( - dt * mu[s] )      # exact solution 

#         # RK2
#         k1        = N[t, s] * mu[s]                    # slope at beginning of interval
#         N_star    = N[t, s] - dt * k1                  # slope at the midpoint of interval
#         k2        = N_star * mu[s]                     # slope at end of interval
#         N[t+1, s] = N[t, s] - (dt * 0.5) * (k1 + k2)   # update value
#         print(N[t+1, s])

# ## ANALYTICAL SOLUTION 
# # initialize analytical solution matrix
# sol = np.zeros([len(time),len(step)])

# # calculate the analytical solution for every age at time t
# for i_t in range(0, len(time)):
#     # Calculate the analytical solution
#     # sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp( - c * time[i_t]) 
#     # sol[i_t,:] = np.exp(-(step - ( time[i_t] + 5))**2) * np.exp( - mu * time[i_t]) 
#     sol[i_t,:] = np.exp(-(step - (0 + 5))**2) * np.exp( - mu * time[i_t]) 


# ## COMPARTISION PLOT BTWN NUMERICAL AND ANALYTICAL
# # get inidices of initial, middle, and last time step
# plot_indices = [0, len(time)  // 2, len(time)  - 1]

# plt.close()
# # plot numerical and analytical solution
# for t_index in plot_indices:
#     plt.plot(step, sol[t_index, :], label=f'Analytical at time {round(time[t_index], 1)  }', linestyle='-')     # analytical 
#     plt.plot(step, N[t_index, :], label=f'Numerical at time {round(time[t_index], 1)  }', linestyle='--')    # numerical 

# # aesthetic of plot
# plt.axhline(y=1, color='r', linestyle='--', label='y=1')
# plt.xlabel('Step')
# plt.ylabel('Population')
# plt.title(f'Population by Step when da = {ds}')
# plt.legend()

# plt.show()