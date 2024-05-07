# Author: Morgan Lavenstein Bendall
# Objective: This calls the function of our model and runs at different ds and dt.

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

# initial conditions
Smax = 15
Tmax = 1
order = 1
c = 0 # constant for mu

Ntest = 5

ds = np.zeros([Ntest]) # order smallest to largest
# ds[4] = 0.1
# ds[3] = 0.2
# ds[2] = 0.3
# ds[1] = 0.4
# ds[0] = 0.5

# # 1 order good - covergence
# ds[0] = 0.06
# ds[1] = 0.07 
# ds[2] = 0.08
# ds[3] = 0.09
# ds[4] = 0.095
# if it too high it might be that the values are too large

# 1 order 
#   - this worked for 2nd order 
# ds[0] = 0.05
# ds[1] = 0.04 
# ds[2] = 0.03
# ds[3] = 0.02
# ds[4] = 0.01

# ds[0] = 0.008
# ds[1] = 0.004 
# ds[2] = 0.002
# ds[3] = 0.001
# ds[4] = 0.0005

# ds[0] = 0.08
# ds[1] = 0.04 
# ds[2] = 0.02
# ds[3] = 0.01
# ds[4] = 0.005

# ds[0] = 0.009
# ds[1] = 0.008 
# ds[2] = 0.007
# ds[3] = 0.006
# ds[4] = 0.005

ds[0] = 0.06
ds[1] = 0.05
ds[2] = 0.04
ds[3] = 0.03
ds[4] = 0.02


# ds[0] = 0.05
# ds[1] = 0.025
# ds[2] = 0.020 
# ds[3] = 0.010
# ds[4] = 0.005


# dt = 0.0001 # da ten times smaller^^

# worked a lot better for order 2 but still not great
# ds[0] = 0.006
# ds[1] = 0.012
# ds[2] = 0.024
# ds[3] = 0.048
# ds[4] = 0.096

cfl = 0.1

dt = np.zeros([Ntest])
dt[0] = cfl * ds[0]
dt[1] = cfl * ds[1]
dt[2] = cfl * ds[2]
dt[3] = cfl * ds[3]
dt[4] = cfl * ds[4]


# dt = np.zeros([ntests])
# dt[0] = 0.006/2
# dt[1] = 0.007/2
# dt[2] = 0.008/2
# dt[3] = 0.009/2
# dt[4] = 0.0095/2

filename = 'ds_convergence/' 


# print('Upwind Order = ' + str(order))
for i in range(len(ds)):

    print('Entering loop ' + str(i))

    # initalize arrays
    Nsize = int(Smax/ds[i]) + 1
    size = np.linspace(0,Smax, Nsize) 


    if isinstance(dt, np.ndarray):
        # initalize arrays
        Ntime = int(Tmax/dt[i]) + 1
        time = np.linspace(0,Tmax, Ntime)

        data = np.zeros([Ntime,Nsize])

        print('CFL: ' + str(round(dt[i]/ds[i], 5)))
        # data = UPW_SPM(size, time, ds[i], dt[i], order, c)
        data = LW_SPM(size, time, ds[i], dt[i], c)
        # LW_SPM(ds[i],dt[i],i,filename)

    else:
        # initalize arrays
        Ntime = int(Tmax/dt) + 1
        time = np.linspace(0,Tmax, Ntime)

        data = np.zeros([Ntime,Nsize])

        print('CFL: ' + str(round(dt/ds[i], 5)))
        # data = UPW_SPM(size, time, ds[i], dt, order, c)
        data = LW_SPM(size, time, ds[i], dt, c)
        # LW_SPM(ds[i],dt,i,filename)

    np.savetxt('ds_convergence/upwind_num_'+ str(i) +'.txt', data)  # save data to file

    print('Loop ' + str(i) + ' Complete.') # Progress update, loop end.


    # analytical solution --------------------------------------------------------------------------------
    sol = np.zeros([len(time),len(size)])
    
    for i_t in range(0, len(time)):
        # Calculate the analytical solution
        sol[i_t,:] = np.exp(-(size - ( time[i_t] + 5))**2) * np.exp( - c * time[i_t]) 


    # Plots the numerical solution at initial, middle, and last time steps
    plot_indices = [0, Ntime // 2, Ntime - 1]  # Indices of initial, middle, and last time steps
    for t_index in plot_indices:
        # plt.plot(size, data[t_index, :], label=f'Numerical at time {round(time[t_index], 1)  }', linestyle='--')
        plt.plot(size, sol[t_index, :], label=f'Analytical at time {round(time[t_index], 1)  }', linestyle='-')
        plt.plot(size, data[t_index, :], label=f'Numerical at time {round(time[t_index], 1)  }', linestyle='--')


    plt.axhline(y=1, color='r', linestyle='--', label='y=1')
    plt.xlabel('Step')
    plt.ylabel('Population')
    plt.title(f'Population by Step when ds = {ds[i] }')
    plt.legend()

    if isinstance(dt, np.ndarray):
        # Save the plot to a file -- labels with da values and dt 
        plt.savefig('ds_plot/varied_dt/plot_mu_' + str(c) + '_ds_' + str(ds[i]) + '_dt_' + str(round(dt[i],5)) + '_order_'+ str(order) +'.png', dpi=300)  
    else:
        # Save the plot to a file -- labels with da values and dt 
        plt.savefig('ds_plot/fixed_dt/plot_mu_' + str(c) + '_ds_' + str(ds[i]) + '_dt_' + str(dt) + '_order_'+ str(order) +'.png', dpi=300) 

    plt.show()



# Plot the convergence xs
if isinstance(dt, np.ndarray):
    Norm = convergence_dt_plt(Smax, Tmax, ds, dt, order, c) 
else:
    Norm = convergence_da_plt(Smax, Tmax, ds, dt, order, c)


# Checks conservation, returns norm and order of conservation
conservation = conservation_plt(Ntest, time, ds, c, Smax, Tmax, dt, order)

# Makes latex table
tabulate_conv(dt, ds, Norm[0], Norm[1], Norm[2], Norm[3], conservation[0], conservation[1])


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