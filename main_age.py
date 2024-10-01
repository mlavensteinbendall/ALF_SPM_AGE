# Author: Morgan Lavenstein Bendall
# Objective: This calls the function of our model and runs at different da and dt.


# Tasks:
#   - need to change file name from upwind to lax-wendroff or something else

import numpy as np
from function_upwind_age    import UPW_SPM
# from convergence_da         import convergence_da_plt
from function_conservation  import conservation_plt
from convergence_dt         import convergence_dt_plt
# from convergence_dt_no_exact_sol         import convergence_dt_plt
from convergence_da         import convergence_da_plt
import matplotlib.pyplot    as plt 
from print_tab_conv         import tabulate_conv
from function_LW            import LW_SPM
# from function_LW_steve      import LW_SPM
import timeit


start = timeit.default_timer()

## INITIAL CONDITIONS
Amax = 15       # max age
Tmax = 0.5      # max time
order = 2       # order of method
c = 0           # constant for mu
Ntest = 5       # number of cases

# folder = 'LW-EX_mu_' + str(c)       # Name of folder for test
folder = 'LW-RK2_mu_' + str(c)       # Name of folder for test

# need to chose da and dt so that the last value in the array are Amax and Tmax
da = np.zeros([Ntest]) # order smallest to largest

# Values that work with: Amax = 15 and Tmax = 0.5 
# da[0] = 0.15
# da[1] = 0.1
# da[2] = 0.05
# da[3] = 0.03
# da[4] = 0.02

# da[0] = 0.15
# da[1] = 0.075
# da[2] = 0.05
# da[3] = 0.025
# da[4] = 0.0125


# da[0] = 0.15
# da[1] = da[0] / 2
# da[2] = da[0] / 4
# da[3] = da[0] / 6
# da[4] = da[0] / 8

# da[0] = 0.15
# da[1] = 0.12
# da[2] = 0.06
# da[3] = 0.04
# da[4] = 0.0375


# da[0] = 0.3
# da[1] = 0.15
# da[2] = 0.1
# da[3] = 0.075
# da[4] = 0.06

# dt = 0.001

# da[0] = dt / 0.1
# da[1] = dt / 0.2
# da[2] = dt / 0.4
# da[3] = dt / 0.5
# da[4] = dt / 0.8


# vary da and dt cases:
da[0] = 0.1
da[1] = 0.05
da[2] = 0.025
da[3] = 0.0125
da[4] = 0.00625

dt = 0.5 * da

# dt = 0.00001

# dt = 0.001



# dt = 0.0001 * da
# dt = 0.02
# dt = 0.00001 # da ten times smaller^^
# dt = 0.0001 # this works for da convergence

filename = 'da_convergence/' 


print('Lax-Wendroff Order = ' + str(order))

# Using the given da and dt values, this loop calculates the numerical solution, solve the analytical 
# solution, and plots the numerical vs. analytical solution. 
# BEGIN LOOP ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for i in range(len(da)):

    print('Entering loop ' + str(i))                # progress update, loop began

    # initalize arrays
    age = np.arange(0, Amax + da[i], da[i])       # array from 0 to Amax
    Nage = len(age)                               # number of elements in age
    print("age:", age[-1])                        # check that last element is Amax

    ## NUMERICAL SOLUTION 
    # IF da and dt are varied, do this -----------------------------------------------------------
    if isinstance(dt, np.ndarray):

        # initalize arrays
        time = np.arange(0,Tmax + dt[i], dt[i])     # array from 0 to Tmax
        Ntime = len(time)                           # number of elements in time
        print("Time:", time[-1])                    # check that last element is Tmax

        # initalize data matrix
        data = np.zeros([Ntime,Nage])

        # print CFL 
        print('CFL: ' + str(round(dt[i]/da[i], 5)))   

        # calculate solution
        # data = UPW_SPM(age, time, da[i], dt[i], order, c)    # upwind method
        data = LW_SPM(age, time, da[i], dt[i], c)              # lax-wendroff method
    
    # ELSE da is varied and dt is constant, do this ------------------------------------------------
    else:
        # initalize arrays
        time = np.arange(0,Tmax + dt, dt)           # array from 0 to Tmax
        Ntime = len(time)                           # number of elements in time
        print("Time:", time[-1])                    # check that the last element is Tmax

        # initialize matrix
        data = np.zeros([Ntime,Nage])

        # print CFL
        print('CFL: ' + str(round(dt/da[i], 5)))

        # calculate solution
        # data = UPW_SPM(age, time, da[i], dt, order, c)       # upwind method
        data = LW_SPM(age, time, da[i], dt, c)                 # lax-wendroff method


    # Save data to a file --------------------------------------------------------------------------
    np.savetxt('da_convergence/num_'+ str(i) +'.txt', data)     # save data to file
    
    print('Loop ' + str(i) + ' Complete.')                      # progress update, loop end


    ## ANALYTICAL SOLUTION 
    # initialize analytical solution matrix
    sol = np.zeros([len(time),len(age)])
    
    # calculate the analytical solution for every age at time t
    for i_t in range(0, len(time)):
        # Calculate the analytical solution
        sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp( - c * time[i_t]) 


    ## COMPARTISION PLOT BTWN NUMERICAL AND ANALYTICAL
    # get inidices of initial, middle, and last time step
    plot_indices = [0, Ntime // 2, Ntime - 1]

    plt.close()
    # plot numerical and analytical solution
    for t_index in plot_indices:
        plt.plot(age, sol[t_index, :], label=f'Analytical at time {round(time[t_index], 1)  }', linestyle='-')     # analytical 
        plt.plot(age, data[t_index, :], label=f'Numerical at time {round(time[t_index], 1)  }', linestyle='--')    # numerical 

    # aesthetic of plot
    plt.axhline(y=1, color='r', linestyle='--', label='y=1')
    plt.xlabel('Step')
    plt.ylabel('Population')
    plt.title(f'Population by Step when da = {da[i] }')
    plt.legend()

    # save plots to folder
    if isinstance(dt, np.ndarray):
        # Save the plot to a file -- labels with da values and dt 
        # plt.savefig('da_plot/varied_dt/lw-ex_plot_mu_' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) +'.png', dpi=300)  
        plt.savefig('da_plot/' + folder + '/varied_dt/lw-ex_plot_mu_' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) + '.png', dpi=300)  
    else:
        # Save the plot to a file -- labels with da values and dt 
        plt.savefig('da_plot/' + folder + '/fixed_dt/lw-ex_plot_mu_' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(dt) + '_order_'+ str(order) + '.png', dpi=300) 

    # show plot
    # plt.show()

    # # error check -- using this with Shilpa's matlab code to make sure we are getting the same values
    # print(data[-1, 99:109])


# END LOOP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pervents corrupting convergence plot
plt.close()

## CONVERGENCE ------------------------------------------------------------------------------------------
# Calculate and plot the convergence, returns an matrix with Norm2, L2norm, NormMax, and LMaxnorm
if isinstance(dt, np.ndarray):
    Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Amax, Tmax, da, dt, order, c, folder) 
else:
    Norm2, L2norm, NormMax, LMaxnorm = convergence_da_plt(Amax, Tmax, da, dt, order, c, folder)

## TOTAL POPULATION ERROR --------------------------------------------------------------------------------
# Checks conservation, returns norm and order of conservation
plt.close()
Norm1, L1norm = conservation_plt(Ntest, time, da, c, Amax, Tmax, dt, order, folder)


## PRINT NORMS --------------------------------------------------------------------------------------------
# print latex table
tabulate_conv(dt, da, Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm, folder, c)

# # print excel compatible table
# if isinstance(dt, np.ndarray):
#     for i in range(len(da)):
#         print(f"{dt[i]}, {da[i]}, {Norm2[i]}, {L2norm[i]}, {NormMax[i]}, {LMaxnorm[i]}, {Norm1[i]}, {L1norm[i]}")
# else:
#     for i in range(len(da)):
#         print(f"{dt}, {da[i]}, {Norm2[i]}, {L2norm[i]}, {NormMax[i]}, {LMaxnorm[i]}, {Norm1[i]}, {L1norm[i]}")


# Print the intial condition --------------------------------------------------------------------------------
# # analytical solution 
# sol = np.zeros([len(age)])
# sol = np.exp(-(age - (5))**2) 

# plt.plot(age, sol)
# plt.xlabel('Age')
# plt.ylabel('Population')
# plt.title(f'Initial Condition by Ages (t = 0)')
# plt.legend()
# plt.show()

stop = timeit.default_timer()

print('Time: ', stop - start)