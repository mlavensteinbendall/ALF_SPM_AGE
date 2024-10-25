# Author: Morgan Lavenstein Bendall
# Objective: This calls the function of our model and runs at different da and dt.


# Tasks:
#   - need to change file name from upwind to lax-wendroff or something else

import numpy as np
from old.function_upwind_age    import UPW_SPM
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
from old.RK2_test             import RK2_function

from function_mortality import mortality


start = timeit.default_timer()

## INITIAL CONDITIONS
Amax = 15       # max age
Tmax = 0.5      # max time
order = 2       # order of method
Ntest = 5       # number of cases


# Mortality set up
m = 1              # constant for mu
c = round(m)
constant = False    # True for constant mu, False for function mu
b = 0              # y-intercept

# Folder for plots and convergence tests
if constant == True:
    # folder = 'LW-EX_mu_' + str(c)       # Name of folder for test
    folder = 'LW-RK2_mu_' + str(c)       # Name of folder for test
else:
    folder = 'LW-RK2_mu_linear_decay_slope_' + str(c)

# LW-RK2_mu_linear_decay

# need to chose da and dt so that the last value in the array are Amax and Tmax
da = np.zeros([Ntest]) # order smallest to largest

# # vary da and dt cases:

da[0] = 0.1
da[1] = 0.1
da[2] = 0.1
da[3] = 0.1
da[4] = 0.1


dt = np.zeros([Ntest]) # order smallest to largest

dt[0] = 0.5 * 0.1
dt[1] = 0.5 * 0.05
dt[2] = 0.5 * 0.025
dt[3] = 0.5 * 0.0125
dt[4] = 0.5 * 0.00625

# dt = 0.5 * da

# dt = 0.00001

# dt = 0.001



# dt = 0.0001 * da
# dt = 0.02
# dt = 0.00001 # da ten times smaller^^
# dt = 0.0001 # this works for da convergence

filename = 'da_convergence/' 


# Using the given da and dt values, this loop calculates the numerical solution, solve the analytical 
# solution, and plots the numerical vs. analytical solution. 
# BEGIN LOOP ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for i in range(len(da)):

    print('Entering loop ' + str(i))                # progress update, loop began

    # initalize arrays
    age = np.arange(0, Amax + da[i], da[i])       # array from 0 to Amax
    Nage = len(age)                               # number of elements in age
    print("age:", age[-1])                        # check that last element is Amax

    mu = np.zeros(Nage)
    mu = mortality(Amax, age, m, b, constant)
    print(mu)

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
        data = RK2_function(age, time, da[i], dt[i], mu)
    
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
        data = RK2_function(age, time, da[i], dt, mu)


    # Save data to a file --------------------------------------------------------------------------
    np.savetxt('da_convergence/num_'+ str(i) +'.txt', data)     # save data to file
    
    print('Loop ' + str(i) + ' Complete.')                      # progress update, loop end


    ## ANALYTICAL SOLUTION 
    # initialize analytical solution matrix
    sol = np.zeros([len(time),len(age)])
    
    # calculate the analytical solution for every age at time t
    for i_t in range(0, len(time)):
        # Calculate the analytical solution
        sol[i_t,:] = np.exp(-(age - 5)**2) * np.exp( - mu * time[i_t])               # no advection



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


    plt.show()

# END LOOP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Pervents corrupting convergence plot
plt.close()

## CONVERGENCE ------------------------------------------------------------------------------------------
# Calculate and plot the convergence, returns an matrix with Norm2, L2norm, NormMax, and LMaxnorm
if isinstance(dt, np.ndarray):
    Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Amax, Tmax, da, dt, order, m, b, constant, folder) 
else:
    Norm2, L2norm, NormMax, LMaxnorm = convergence_da_plt(Amax, Tmax, da, dt, order, m, b, constant, folder)

## TOTAL POPULATION ERROR --------------------------------------------------------------------------------
# Checks conservation, returns norm and order of conservation
plt.close()
Norm1, L1norm = conservation_plt(Ntest, time, da, c, Amax, Tmax, dt, order, folder)   # only works for constant 


## PRINT NORMS --------------------------------------------------------------------------------------------
# print latex table
tabulate_conv(dt, da, Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm, folder, c)


stop = timeit.default_timer()

print('Time: ', stop - start)