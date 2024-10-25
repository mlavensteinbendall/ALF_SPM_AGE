import numpy as np # Numpy for numpy
import math
import matplotlib.pyplot as plt
from tabulate import tabulate

def trapezoidal_rule(fx, dx):
    """Performs trapezoidal rule
    
    Args:
        fx  (array):    A list of the population at different steps.
        dx  (int):      The partition of steps.
        
    Returns:
        result  (array): Represents the time
    """

    fx_sum = np.sum(fx[1:-1])

    result = dx * ( (fx[0] + fx[-1]) / 2 + fx_sum)

    return result



def conservation_plt(Ntest, time, ds, c, Smax, Tmax, dt, order, folder):
    """Performs trapezoidal rule
    
    Args:
        Ntest  (array):    A list of the population at different steps.
        time    ():
        ds      (int):
        c       (int):
        Smax    ():
        Tmax    ():
        dt      ():
        order   ():
        
    Returns:
        Norm1  (array):     A list of 1-norm errors
        L1norm  (array):    A list of 1-norm orders.
    """

    totalPop_num = np.zeros([5]) 

    # totalPop_sol = 0 # for advection with mu = 0 it will remain the same everywhere

    Norm1 = np.zeros([5])
    L1norm = np.zeros([5])

    print(time[-1])

    # Exact solution 
    if c == 0:
        totalPop_sol = -0.5 * np.pi**(0.5) * ( math.erf(5 - Smax) - math.erf(5) ) # solution
    else:
        totalPop_sol = -0.5 * np.pi**(0.5) * np.exp(-c * Tmax) * ( math.erf(Tmax + 5 - Smax) - math.erf(Tmax + 5) ) # solution

        print('Exact total pop = ' + str(totalPop_sol))

    # Calculate the total population using trapezoidal rule
    for i in range(5):

        data = np.loadtxt('da_convergence/num_' + str(i) + '.txt') # Load in relevant data.
        # print(data)

        totalPop_num[i] = trapezoidal_rule( data[-1,:],     ds[i])

        # print( 'ds = ' + str(ds[i]) + ' & dt = ' + str(dt[i]) + ' : total pop = ' + str(totalPop_num[i]))
        # np.set_printoptions(precision=15)
        # print(data[-1,119:125])
    

        # Nstep = int(Smax/ds[i]) + 1   # total number of steps

        Norm1[i]    = np.abs( totalPop_num[i] - totalPop_sol )

    for ii in range(0, Ntest - 1):
        L1norm[ii+1] = np.log( Norm1[ii]   / Norm1[ii+1] )   / np.log( ds[ii] / ds[ii+1] )


    for i in range(0, 5):

        print('For ds ='    + str( round( ds[i],10      ) ) )
        print('Norm1 (abs error): '   + str( round( Norm1[i], 20  ) ) )
        if i > 0:
            print('L1 q order: ' + str( round( L1norm[i], 10  ) ) )


    # Plot absolute error oftotal population over time
    # plt.figure(figsize=(8, 6))  # Adjust the width and height as needed
    plt.loglog(ds, Norm1, label='Norm1')
    plt.loglog(ds, ds**(order), label=f'order-{(order) }')
    plt.xlabel('ds')
    plt.ylabel('Absolute Error')
    plt.title('Error of Total Population')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, np.round(ds, 3) ))

    # # Save the plot to a file -- labels with da values and dt 
    # if isinstance(dt, np.ndarray):
    #     dt_values_str = '_'.join(map(str, np.round(dt, 3)))

    #     plt.savefig('da_plot/'+ folder +'/varied_dt/lw-ex_plot_totPop_mu_' + str(c) + '_ds_' + ds_values_str + '_dt_' + dt_values_str + '.png', dpi=300)  

    # else:
    #     plt.savefig('da_plot/'+ folder +'/fixed_dt/lw-ex_plot_totPop_mu_' + str(c) + '_ds_' + ds_values_str + '_dt_' + str(dt) + '.png', dpi=300)  

    plt.show()


    # combine = [Norm1, L1norm]

    return Norm1, L1norm



    # # Plot total population over time
    # plt.plot(ds, totalPop_num)
    # plt.xlabel('ds')
    # plt.ylabel('total pop')
    # plt.title('Total Population at final time')
    # plt.show()