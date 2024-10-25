# Author: Morgan Lavenstein Bendall
# Objective: This calls the function of our model and runs at different da and dt.

import numpy as np
from old.function_upwind_age    import UPW_SPM
# from convergence_da         import convergence_da_plt
from function_conservation  import conservation_plt
# from convergence_dt         import convergence_dt_plt
from convergence_dt_no_exact_sol         import convergence_dt_plt
from convergence_da         import convergence_da_plt
import matplotlib.pyplot    as plt 
from print_tab_conv         import tabulate_conv
from function_LW            import LW_SPM
import timeit
from old.RK2_test             import RK2_function

from function_mortality import mortality


start = timeit.default_timer()

## INITIAL CONDITIONS
Amax = 30       # max age
Tmax = 5     # max time
order = 2       # order of method
Ntest = 5       # number of cases


# Mortality set up
m = 0.5            # constant for mux
c = round(m)
constant = True    # True for constant mu, False for function mu
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
da[1] = 0.05
da[2] = 0.025
da[3] = 0.0125
da[4] = 0.00625

# da[0] = 0.001
# da[1] = 0.0005
# da[2] = 0.00025
# da[3] = 0.000125
# da[4] = 0.0000625


dt = np.zeros([Ntest]) # order smallest to largest

# dt[0] = 0.5 * 0.1
# dt[1] = 0.5 * 0.05
# dt[2] = 0.5 * 0.025
# dt[3] = 0.5 * 0.0125
# dt[4] = 0.5 * 0.00625

# dt = 1 * da
dt = 0.5 * da
# dt = da / 8 


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
    # print(mu)

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
        data = LW_SPM(age, time, da[i], dt[i], mu)              # lax-wendroff method
    
    # ELSE da is varied and dt is constant, do this ------------------------------------------------
    else:
        # initalize arrays
        time = np.arange(0, Tmax + dt, dt)           # array from 0 to Tmax
        Ntime = len(time)                           # number of elements in time
        print("Time:", time[-1])                    # check that the last element is Tmax

        # initialize matrix
        data = np.zeros([Ntime, Nage])

        # print CFL
        print('CFL: ' + str(round(dt/da[i], 5)))

        # calculate solution
        data = LW_SPM(age, time, da[i], dt, mu)                 # lax-wendroff method


    # Save data to a file --------------------------------------------------------------------------
    np.savetxt('da_convergence/num_'+ str(i) +'.txt', data)     # save data to file
    
    print('Loop ' + str(i) + ' Complete.')                      # progress update, loop end


    # ## PLOT MORTALITY RATe
    # plt.plot(age, mu)
    # plt.xlabel('age')
    # plt.ylabel('Mortality rate')
    # plt.title(f'Mortality rate for age')
    # plt.show()

    # plt.close()


    # ## ANALYTICAL SOLUTION 
    # # initialize analytical solution matrix
    # sol = np.zeros([len(time),len(age)])
    
    # # calculate the analytical solution for every age at time t
    # for i_t in range(0, len(time)):
    #     # Calculate the analytical solution
    #     # sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp( - mu * time[i_t])     # with advection -- CONSTANT
    #     sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- m * (age )* time[i_t] + 0.5 * m * (time[i_t])**2)     # with advection -- NON CONSTANT

    #     # sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- (time[i_t] + 200 * np.log(age**2 + 400) - 200 * np.log((age - time[i_t])**2 +400))) # with advection -- hill function

    #     # sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- (20 * np.log(age**2 + 400) - 20 * np.log((age - time[i_t])**2 +400))) # with advection -- hill function



    # # COMPARTISION PLOT BTWN NUMERICAL AND ANALYTICAL
    # # get inidices of initial, middle, and last time step
    # plot_indices = [0, Ntime // 2, Ntime - 1]

    # plt.close()
    # # plot numerical and analytical solution
    # for t_index in plot_indices:
    #     plt.plot(age, sol [t_index, :], label=f'Analytical at time {round(time[t_index], 1)  }', linestyle='-')     # analytical 
    #     plt.plot(age, data[t_index, :], label=f'Numerical at time  {round(time[t_index], 1)  }', linestyle='--')    # numerical 

    # # aesthetic of plot
    # plt.axhline(y=1, color='r', linestyle='--', label='y=1')
    # plt.xlabel('Step')
    # plt.ylabel('Population')
    # plt.title(f'Population by Step when da = {da[i] }')
    # plt.legend()

    # # save plots to folder
    # if isinstance(dt, np.ndarray):
    #     # Save the plot to a file -- labels with da values and dt 
    #     # plt.savefig('da_plot/varied_dt/lw-ex_plot_mu_' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) +'.png', dpi=300)  
    #     if constant == True:
    #         plt.savefig('da_plot/' + folder + '/varied_dt/lw-RK2_plot_mu_' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) + '.png', dpi=300)  
    #     else:
    #          plt.savefig('da_plot/' + folder + '/varied_dt/lw-RK2_plot_mu_slope' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) + '.png', dpi=300)  
    # else:
    #     # Save the plot to a file -- labels with da values and dt 
    #     if constant == True:
    #         plt.savefig('da_plot/' + folder + '/fixed_dt/lw-RK2_plot_mu_' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(dt) + '_order_'+ str(order) + '.png', dpi=300) 
    #     else:
    #          plt.savefig('da_plot/' + folder + '/fixed_dt/lw-RK2_plot_mu_slope' + str(c) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) + '.png', dpi=300)  
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
    # Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Amax, Tmax, da, dt, order, m, b, constant, folder) 
    Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Amax, da, dt, order)
else:
    Norm2, L2norm, NormMax, LMaxnorm = convergence_da_plt(Amax, Tmax, da, dt, order, m, b, constant, folder)


## TOTAL POPULATION ERROR --------------------------------------------------------------------------------
# Checks conservation, returns norm and order of conservation
plt.close()
Norm1, L1norm = conservation_plt(Ntest, time, da, c, Amax, Tmax, dt, order, folder)   # only works for constant 


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