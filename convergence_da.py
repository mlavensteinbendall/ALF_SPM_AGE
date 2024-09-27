import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def convergence_da_plt(age_max, time_max, da, dt, order, c, folder):
    """Calculates the convergence for varying ds and constant dt.
    
    Args:
        age_max (int):      The max number of age-steps.
        time_max(float):    The max number of time-steps.
        da      (float):    The partitions for age-steps.
        dt      (float):    The partitions for time-steps.
        order   (int):      The order of the method.
        c       (int):      The constant of death rate.
        
    Returns:
        Norm2   (array):    A list of the 2-norms.
        L2norm  (array):    A list of the order of the 2-norms.
        NormMax (array):    A list of the infinity-norms.
        LMaxnorm(array):    A list of the order of the infinity-norms.
    """

    n = int(time_max/dt) + 1 # Time-step of comparison.
    Tend = n*dt # Get the associated timepoint value.

    # Tend = time_max
    Ntest = len(da)
    print(Ntest)

    Norm2 = np.zeros([Ntest])
    NormMax = np.zeros([Ntest])

    L2norm = np.zeros([Ntest])
    LMaxnorm = np.zeros([Ntest])


    for i in range(0, Ntest):

        # initialize values
        age = np.arange(0, age_max + da[i], da[i])      # array from 0 to age_max
        Nage = len(age)                                 # number of elements in age

        data = np.zeros([int(time_max/da[i]), Nage])    # initialize matrix for numerical solution
        sol = np.zeros([Nage])                          # initialize array for analytical solution

        # Numerical solution -- download relevent data
        data = np.loadtxt('da_convergence/num_' + str(i) + '.txt') # Load in relevant data.

        # Analyticial solution -- changes for what ds is
        sol = np.exp(-(age - ( Tend + 5))**2) * np.exp(-c * Tend)        

        
        # # plt data vs sol
        # plt.plot(step,data[-1,:]) 
        # plt.plot(step,sol) # looks right
        # plt.show()

        # Calculate the norm 
        Norm2[i]    = ( ( 1 / Nage ) * np.sum( np.abs( data[-1,:] - sol[:] ) **2 ) ) **0.5  # L2 error.
        NormMax[i]  = np.max( np.abs( data[-1,:] - sol[:] ) )                               # L-Max error.


    # Iterate to calculates the L norms -- comparing with the last (Note: ds are decressing)
    for ii in range(0, Ntest - 1):
        L2norm[ii+1]    = np.log( Norm2[ii+1]   / Norm2[ii] )   / np.log( da[ii+1] / da[ii] )
        LMaxnorm[ii+1]  = np.log( NormMax[ii+1] / NormMax[ii] ) / np.log( da[ii+1] / da[ii] )


    # Print error and order for each combination of ds and dt
    for i in range(0, Ntest):

        print('For ds ='            + str( round( da[i],        10  ) ) )
        print('Norm 2 error: '      + str( round( Norm2[i],     10  ) ) )
        print('Norm inf error: '    + str( round( NormMax[i],   10  ) ) )

        if i > 0:
            print('L2 q order: '    + str( round( L2norm[i-1]   , 10    ) ) ) # L2 q estimate.
            print('LMax q order: '  + str( round( LMaxnorm[i-1] , 10    ) ) ) # L-Max q estimate.
            print(' ')


    # Plot the log-log for the errors.
    plt.loglog(da, Norm2, label='Norm2')
    plt.loglog(da, NormMax, label='NormMax')
    plt.loglog(da, da**(order), label=f'order-{order }')

    plt.xlabel(r'$\Delta a$')
    plt.ylabel('Norm')
    plt.title('Convergence based on ' + r'$\Delta a$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, da))

    # Save the plot to a file -- labels with da values and dt 
    plt.savefig('da_plot/' + folder + '/fixed_dt/lw-ex_plot_conv_mu_'+ str(c) + '_ds_'+ ds_values_str + f'_dt_{dt }' + '.png', dpi=300)  

    plt.show()

    return Norm2, L2norm, NormMax, LMaxnorm

    # # Plot the Analytical and Numerical Solution
    # plt.plot(size, sol, label='Analytical', linestyle='solid')
    # plt.plot(size, data[-1,:], label='Numerical', linestyle='--')
    # # plt.axvline(x=0.4, color='r', linestyle='--', label='x=1')
    # plt.xlabel('Size')
    # plt.ylabel('Population')
    # plt.title('Population Based on Size at time = 0')
    # plt.legend()
    # # plt.ylim(-.2, 1.4)  # Set y-axis limits from 0 to 12
    # # plt.savefig('plots/pop_plot-time' + str(n) + '-ds_'+ str(ds_index) +'.png') # Save the plot
    # plt.show()