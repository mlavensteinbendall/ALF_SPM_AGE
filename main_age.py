# Author: Morgan Lavenstein Bendall
# Objective: This calls the function of our model and runs at different ds and dt.


# Tasks:
#   - need to change file name from upwind to lax-wendroff or something else

import numpy as np
from function_upwind_age    import UPW_SPM
# from convergence_ds         import convergence_ds_plt
from function_conservation  import conservation_plt
from convergence_dt         import convergence_dt_plt
# from convergence_dt_no_exact_sol         import convergence_dt_plt
from convergence_da         import convergence_da_plt
import matplotlib.pyplot    as plt 
from print_tab_conv         import tabulate_conv
from function_LW            import LW_SPM
# from function_LW_steve      import LW_SPM

## INITIAL CONDITIONS
Smax = 15       # max age
Tmax = 0.5      # max time
order = 2       # order of method
c = 1           # constant for mu
Ntest = 5       # number of cases

# need to chose ds and dt so that the last value in the array are Smax and Tmax
ds = np.zeros([Ntest]) # order smallest to largest

# Values that work with: Smax = 15 and Tmax = 0.5 
# ds[0] = 0.15
# ds[1] = 0.1
# ds[2] = 0.05
# ds[3] = 0.03
# ds[4] = 0.02

# ds[0] = 0.15
# ds[1] = 0.075
# ds[2] = 0.05
# ds[3] = 0.025
# ds[4] = 0.0125

# ds[0] = 0.3
# ds[1] = 0.15
# ds[2] = 0.1
# ds[3] = 0.075
# ds[4] = 0.06

# ds[0] = 0.15
# ds[1] = 0.12
# ds[2] = 0.06
# ds[3] = 0.04
# ds[4] = 0.0375

# dt = 0.001

# vary ds and dt cases:
ds[0] = 0.0625
ds[1] = 0.05
ds[2] = 0.025
ds[3] = 0.0125
ds[4] = 0.003125

dt = 0.5 * ds

# dt = 0.0001 * ds
# dt = 0.02
# dt = 0.00001 # da ten times smaller^^
# dt = 0.0001 # this works for da convergence

filename = 'ds_convergence/' 


print('Lax-Wendroff Order = ' + str(order))

# Using the given ds and dt values, this loop calculates the numerical solution, solve the analytical 
# solution, and plots the numerical vs. analytical solution. 
# BEGIN LOOP ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for i in range(len(ds)):

    print('Entering loop ' + str(i))                # progress update, loop began

    # initalize arrays
    size = np.arange(0, Smax + ds[i], ds[i])        # array from 0 to Smax
    Nsize = len(size)                               # number of elements in size
    print("Size:", size[-1])                        # check that last element is Smax

    ## NUMERICAL SOLUTION 
    # IF ds and dt are varied, do this -----------------------------------------------------------
    if isinstance(dt, np.ndarray):

        # initalize arrays
        time = np.arange(0,Tmax + dt[i], dt[i])     # array from 0 to Tmax
        Ntime = len(time)                           # number of elements in time
        print("Time:", time[-1])                    # check that last element is Tmax

        # initalize data matrix
        data = np.zeros([Ntime,Nsize])

        # print CFL 
        print('CFL: ' + str(round(dt[i]/ds[i], 5)))   

        # calculate solution
        # data = UPW_SPM(size, time, ds[i], dt[i], order, c)    # upwind method
        data = LW_SPM(size, time, ds[i], dt[i], c)              # lax-wendroff method
    
    # ELSE ds is varied and dt is constant, do this ------------------------------------------------
    else:
        # initalize arrays
        time = np.arange(0,Tmax + dt, dt)           # array from 0 to Tmax
        Ntime = len(time)                           # number of elements in time
        print("Time:", time[-1])                    # check that the last element is Tmax

        # initialize matrix
        data = np.zeros([Ntime,Nsize])

        # print CFL
        print('CFL: ' + str(round(dt/ds[i], 5)))

        # calculate solution
        # data = UPW_SPM(size, time, ds[i], dt, order, c)       # upwind method
        data = LW_SPM(size, time, ds[i], dt, c)                 # lax-wendroff method


    # Save data to a file --------------------------------------------------------------------------
    np.savetxt('ds_convergence/num_'+ str(i) +'.txt', data)     # save data to file
    
    print('Loop ' + str(i) + ' Complete.')                      # progress update, loop end


    ## ANALYTICAL SOLUTION 
    # initialize analytical solution matrix
    sol = np.zeros([len(time),len(size)])
    
    # calculate the analytical solution for every size at time t
    for i_t in range(0, len(time)):
        # Calculate the analytical solution
        sol[i_t,:] = np.exp(-(size - ( time[i_t] + 5))**2) * np.exp( - c * time[i_t]) 


    ## COMPARTISION PLOT BTWN NUMERICAL AND ANALYTICAL
    # get inidices of initial, middle, and last time step
    plot_indices = [0, Ntime // 2, Ntime - 1]

    # plot numerical and analytical solution
    for t_index in plot_indices:
        plt.plot(size, sol[t_index, :], label=f'Analytical at time {round(time[t_index], 1)  }', linestyle='-')     # analytical 
        plt.plot(size, data[t_index, :], label=f'Numerical at time {round(time[t_index], 1)  }', linestyle='--')    # numerical 

    # aesthetic of plot
    plt.axhline(y=1, color='r', linestyle='--', label='y=1')
    plt.xlabel('Step')
    plt.ylabel('Population')
    plt.title(f'Population by Step when ds = {ds[i] }')
    plt.legend()

    # save plots to folder
    if isinstance(dt, np.ndarray):
        # Save the plot to a file -- labels with da values and dt 
        plt.savefig('ds_plot/varied_dt/plot_mu_' + str(c) + '_ds_' + str(ds[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) +'.png', dpi=300)  
    else:
        # Save the plot to a file -- labels with da values and dt 
        plt.savefig('ds_plot/fixed_dt/plot_mu_' + str(c) + '_ds_' + str(ds[i]) + '_dt_' + str(dt) + '_order_'+ str(order) +'.png', dpi=300) 

    # show plot
    plt.show()

    # # error check -- using this with Shilpa's matlab code to make sure we are getting the same values
    # print(data[-1, 99:109])


# END LOOP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## CONVERGENCE ------------------------------------------------------------------------------------------
# Calculate and plot the convergence, returns an matrix with Norm2, L2norm, NormMax, and LMaxnorm
if isinstance(dt, np.ndarray):
    Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Smax, Tmax, ds, dt, order, c) 
else:
    Norm2, L2norm, NormMax, LMaxnorm = convergence_da_plt(Smax, Tmax, ds, dt, order, c)

## TOTAL POPULATION ERROR --------------------------------------------------------------------------------
# Checks conservation, returns norm and order of conservation
Norm1, L1norm = conservation_plt(Ntest, time, ds, c, Smax, Tmax, dt, order)


## PRINT NORMS --------------------------------------------------------------------------------------------
# print latex table
tabulate_conv(dt, ds, Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm)

# # print excel compatible table
# if isinstance(dt, np.ndarray):
#     for i in range(len(ds)):
#         print(f"{dt[i]}, {ds[i]}, {Norm2[i]}, {L2norm[i]}, {NormMax[i]}, {LMaxnorm[i]}, {Norm1[i]}, {L1norm[i]}")
# else:
#     for i in range(len(ds)):
#         print(f"{dt}, {ds[i]}, {Norm2[i]}, {L2norm[i]}, {NormMax[i]}, {LMaxnorm[i]}, {Norm1[i]}, {L1norm[i]}")


# Print the intial condition --------------------------------------------------------------------------------
# # analytical solution 
# sol = np.zeros([len(size)])
# sol = np.exp(-(size - (5))**2) 

# plt.plot(size, sol)
# plt.xlabel('Age')
# plt.ylabel('Population')
# plt.title(f'Initial Condition by Ages (t = 0)')
# plt.legend()
# plt.show()