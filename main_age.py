# Author: Morgan Lavenstein Bendall
# Objective: This calls the function of our model and runs at different da and dt.

import numpy                        as np
import matplotlib.pyplot            as plt 
import timeit
from old.function_upwind_age        import UPW_SPM
# from convergence_da         import convergence_da_plt
# from function_conservation  import conservation_plt
from function_conservation_no_exact_sol import conservation_plt
# from convergence_dt         import convergence_dt_plt
from convergence_dt_no_exact_sol    import convergence_dt_plt
from convergence_da_no_exact_sol    import convergence_da_plt
# from convergence_da         import convergence_da_plt
from print_tab_conv                 import tabulate_conv
from function_LW                    import LW_SPM
from old.RK2_test                   import RK2_function
from function_reproduction          import reproduction

from function_mortality import mortality


start = timeit.default_timer()

## INITIAL CONDITIONS
Amax = 30       # max age
Tmax = 5     # max time
order = 2       # order of method
Ntest = 5       # number of cases


# Mortality set up
m = 0       #1/30            # constant for mux
b = 0              # y-intercept
constant_mortality = True   # True for constant mu, False for function mu
analytical_sol = False
hill_func_mortality = False
linear_slope_func_mortality = False

# Reproduction set up
rep = 1
constant_reproduction = True
linear_reproduction = False

# testing_folder = 'mortality'
testing_folder = 'reproduction'

if constant_mortality == True:
    if m == 0 :
        function_folder = "no_mortality"

    else:
        function_folder = "constant_mortality"

elif hill_func_mortality == True:
    function_folder = "hill_mortality"

elif linear_slope_func_mortality == True:
    function_folder = "linear_mortality"

else:
    function_folder = "gompertz_mortality"

if testing_folder == 'reproduction':

    if constant_reproduction == True:
        if rep == 0:
            function_folder = function_folder + '/' + 'no_reproduction'

        else:
            function_folder = function_folder + '/' + 'constant_reproduction/rep_' + str(rep)

    elif linear_reproduction == True:
        function_folder = function_folder + '/' + 'linear_reproduction'




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

# da[0] = 1
# da[1] = 0.5
# da[2] = 0.25
# da[3] = 0.125
# da[4] = 0.0625


dt = np.zeros([Ntest]) # order smallest to largest

# dt[0] = 0.5 * 0.1
# dt[1] = 0.5 * 0.05
# dt[2] = 0.5 * 0.025
# dt[3] = 0.5 * 0.0125
# dt[4] = 0.5 * 0.00625

# dt = 0.5 * da
# dt = 0.01
# dt = 0.001
# dt = 0.0001
dt = 0.000001
# da = 0.5 * da
# dt = 0.5 * dt

if isinstance(dt, np.ndarray):
    convergence_folder = 'varied_dt'

else:
    convergence_folder = 'fixed_dt'


folder = 'convergence/' + testing_folder + '/' + function_folder + '/' + convergence_folder


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
    mu = mortality(Amax, age, m, b, constant_mortality, linear_slope_func_mortality, hill_func_mortality)

    if i == 0:
        plt.plot(age, mu)
        plt.xlabel('Age')
        plt.ylabel('Mortality Rate')
        plt.title('Age-Specific Mortality Rate')
        if isinstance(dt, np.ndarray): 
            plt.savefig(folder + '/plots/mortality_plot.png', dpi=300)
        else:
            plt.savefig(folder + '/plots/dt_' + str(dt) + '/mortality_plot.png', dpi=300)
        # plt.show()
        plt.close()


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
        data = LW_SPM(age, time, da[i], dt[i], mu, rep, constant_reproduction)              # lax-wendroff method
    
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
        data = LW_SPM(age, time, da[i], dt, mu, rep, constant_reproduction)                 # lax-wendroff method


    # Save data to a file --------------------------------------------------------------------------
    if isinstance(dt, np.ndarray):   
        np.savetxt(folder + '/solutions/num_'+ str(i) +'.txt', data)     # save data to file 

    else:
        np.savetxt(folder + '/solutions/dt_' + str(dt) + '/num_'+ str(i) +'.txt', data)     # save data to file 

    print('Loop ' + str(i) + ' Complete.')                      # progress update, loop end


    
    # calculate the analytical solution for every age at time t
    for i_t in range(0, len(time)):
        # Calculate the analytical solution

        # for i_a in range(0, len(age)):

            # if constant == True:
            #     sol[i_t, i_a] = np.exp(-(age[i_a] - ( time[i_t] + 5))**2) * np.exp( - mu[i_a] * time[i_t])     # with advection -- CONSTANT

            # elif hill_func == True:
            #     sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- (30 * np.log(age**2 + 30**2) - 30 * np.log((age - time[i_t])**2 +30**2))) # with advection -- hill function

            # elif linear_slope_func == True:
            #     sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- m * (age )* time[i_t] + 0.5 * m * (time[i_t])**2)     # with advection -- NON CONSTANT linear slope

        # if i_t > 0:
            # print('type = ' + str(type(sol[i_t, :])))
            # print('shape = ' + str((sol[i_t, :]).shape))
            # print(type(sol[i_t, 0]))
            # print('beepbeep' + str(type(reproduction(sol[i_t, :], 0, da) )))
            # temp = sol[i_t, :]
            # print(temp.shape)
            # sol[, 0] = reproduction(temp, 0, da) 
            # print(type(sol[i_t,0]))
            # print("type of boundry:" + str(type(reproduction(temp, 0, da) )))
        
        
        ## ANALYTICAL SOLUTION 
        if analytical_sol == True:

             # initialize analytical solution matrix
            sol = np.zeros([len(time),len(age)])

            if constant_mortality == True:
                sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp( - mu * time[i_t])     # with advection -- CONSTANT

            elif hill_func_mortality == True:
                sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- (30 * np.log(age**2 + 30**2) - 30 * np.log((age - time[i_t])**2 +30**2))) # with advection -- hill function

            elif linear_slope_func_mortality == True:
                sol[i_t,:] = np.exp(-(age - ( time[i_t] + 5))**2) * np.exp(- m * (age )* time[i_t] + 0.5 * m * (time[i_t])**2)     # with advection -- NON CONSTANT linear slope


    # COMPARTISION PLOT BTWN NUMERICAL AND ANALYTICAL
    # get inidices of initial, middle, and last time step
    plot_indices = [0, Ntime // 2, Ntime - 1]

    plt.close()
    # plot numerical and analytical solution
    for t_index in plot_indices:
        if analytical_sol == True: # check if there are non-zero elements in the array (if it's all zeros, then there was no analytical solution to solve)
            plt.plot(age, sol [t_index, :], label=f'Analytical at time {round(time[t_index], 1)  }', linestyle='-')     # analytical 
        plt.plot(age, data[t_index, :], label=f'Numerical at time  {round(time[t_index], 1)  }', linestyle='--')    # numerical 

    plt.axhline(y=1, color='r', linestyle='--', label='y=1')
    plt.xlabel('Age')
    plt.ylabel('Population')
    if isinstance(dt, np.ndarray):
        # plt.title(f'Population by Step when $\Delta a$ = {da[i] } and $\Delta t$ = {dt[i] }')
        plt.title('Age Distribution of Population (' + r'$\Delta a$' + ' = ' + str(da[i]) + ', ' + r'$\Delta t$' + ' = ' + str(dt[i]) + ')')
    else:
        plt.title('Age Distribution of Population (' + r'$\Delta a$' + ' = ' + str(da[i]) + ', ' + r'$\Delta t$' + ' = ' + str(dt) + ')')
        # plt.title(f'Population by Step when' + r'$\Delta a$' + f' = {da[i] } and ' + r'$\Delta t' + f' = {dt }')
    plt.legend()

    # save plots to folder
    if isinstance(dt, np.ndarray):

        plt.savefig(folder + '/plots/num_' + str(i) + '_da_' + str(da[i]) + '_dt_' + str(round(dt[i],5)) + '.png', dpi=300)
    else:

        plt.savefig(folder + '/plots/dt_' + str(dt) + '/num_' + str(i) + '_da_' + str(da[i]) + '_dt_' + str(dt) + '.png', dpi=300)
    
    # plt.show()        # show plot

    # # error check -- using this with Shilpa's matlab code to make sure we are getting the same values
    # print(data[-1, 99:109])


# END LOOP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Pervents corrupting convergence plot
plt.close()

## CONVERGENCE ------------------------------------------------------------------------------------------
# Calculate and plot the convergence, returns an matrix with Norm2, L2norm, NormMax, and LMaxnorm
if isinstance(dt, np.ndarray):
    # Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Amax, Tmax, da, dt, order, m, b, constant, folder) 
    Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Tmax, da, dt, order, folder)
else:
    # Norm2, L2norm, NormMax, LMaxnorm = convergence_da_plt(Amax, Tmax, da, dt, order, m, b, constant, folder)
    Norm2, L2norm, NormMax, LMaxnorm = convergence_da_plt(Tmax, da, dt, order, folder)


## TOTAL POPULATION ERROR --------------------------------------------------------------------------------
# Checks conservation, returns norm and order of conservation
plt.close()
# Norm1, L1norm = conservation_plt(Ntest, da, m, Amax, Tmax, dt, order, folder, constant, hill_func)   # only works for constant 
Norm1, L1norm = conservation_plt(da, dt, order, folder)


## PRINT NORMS --------------------------------------------------------------------------------------------
# print latex table
tabulate_conv(dt, da, Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm, folder)

# # print excel compatible table
# if isinstance(dt, np.ndarray):
#     for i in range(len(da)):
#         print(f"{dt[i]}, {da[i]}, {Norm2[i]}, {L2norm[i]}, {NormMax[i]}, {LMaxnorm[i]}, {Norm1[i]}, {L1norm[i]}")
# else:
#     for i in range(len(da)):
#         print(f"{dt}, {da[i]}, {Norm2[i]}, {L2norm[i]}, {NormMax[i]}, {LMaxnorm[i]}, {Norm1[i]}, {L1norm[i]}")


stop = timeit.default_timer()

print('Time: ', stop - start)