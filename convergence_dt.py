import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def convergence_dt_plt(Smax, Tmax, ds, dt, order, c, folder):
    """Calculates the convergence for varying ds and dt.
    
    Args:
        Smax    (int):      The max number of steps.
        Tmax    (float):    The max number of time-steps.
        ds      (float):    The partitions for steps.
        dt      (float):    The partitions for time-steps.
        order   (int):      The order of the method.
        c       (int):      The constant of death rate.
        
    Returns:
        Norm2   (array):    A list of the 2-norms.
        L2norm  (array):    A list of the order of the 2-norms.
        NormMax (array):    A list of the infinity-norms.
        LMaxnorm(array):    A list of the order of the infinity-norms.
    """

    Ntest = len(ds)                 # number of cases

    # initialize the arrays for the norms and norm orders to zero
    Norm2 = np.zeros([Ntest])       # 2 norm
    NormMax = np.zeros([Ntest])     # infinity norm
    L2norm = np.zeros([Ntest])      # order for 2 norm
    LMaxnorm = np.zeros([Ntest])    # order for infinity norm


    # Iterate over the number of tests (0 to Ntest)
    for i in range(0, Ntest):

        # initialize values
        step = np.arange(0, Smax + ds[i], ds[i])    # array from o to Smax
        Nstep = len(step)                           # number of elements in step

        Ntime = int(Tmax/dt[i])     # Time-step of comparison.
        Tend = Ntime * dt[i]        # Get the associated timepoint value.

        # data = np.zeros([int(Tmax/ds[i]), Nstep])       # initialize matrix for numerical solution
        data = np.zeros([Ntime, Nstep])                   # initialize matrix for numerical solution
        sol = np.zeros([Nstep])                           # initialize array for analytical solution


        # Numerical solution -- download relevent data
        data = np.loadtxt('da_convergence/num_' + str(i) + '.txt') 

        # Analyticial solution -- changes for what ds is
        sol = np.exp(-(step - ( Tend + 5))**2) * np.exp(-c * Tend)         
        
        # # plt data vs sol
        # plt.plot(step,data[-1,:]) 
        # plt.plot(step,sol) # looks right
        # plt.show()

        # Calculate the norm 
        Norm2[i]    = ( ( 1 / Nstep ) * np.sum( np.abs( data[-1,:] - sol[:] ) **2 ) ) **0.5     # L2 error
        NormMax[i]  = np.max( np.abs( data[-1,:] - sol[:] ) )                                   # Lmax error


    # Iterate to calculates the L norms -- comparing with the last (Note: ds and dt are decressing with same CFL value)
    for ii in range(0, Ntest - 1):
        L2norm[ii+1]    = np.log( Norm2[ii+1]   / Norm2[ii] )   / np.log( dt[ii+1] / dt[ii] )   # order from L2
        LMaxnorm[ii+1]  = np.log( NormMax[ii+1] / NormMax[ii] ) / np.log( dt[ii+1] / dt[ii] )   # order from Lmax


    # Print error and order for each combination of ds and dt
    for i in range(0, Ntest):

        print('For ds ='    + str( round( ds[i],10      ) ) + ' and dt ='    + str( round( dt[i],10      ) ) )
        print('Norm 2 error: '   + str( round( Norm2[i], 10  ) ) )
        print('Norm inf error: ' + str( round( NormMax[i], 10) ) )

        if i > 0:
            print('L2 q order: '     + str( round( L2norm[i-1]   , 10    ))) # L2 q estimate.
            print('LMax q order: '   + str( round( LMaxnorm[i-1] , 10    ))) # L-Max q estimate.
            print(' ')


    plt.clf()
    # Plot the log-log for the errors.
    plt.loglog(ds, Norm2, label='Norm2')
    plt.loglog(ds, NormMax, label='NormMax')
    plt.loglog(ds, ds**(order), label=f'order-{order }')

    plt.xlabel(r'$\Delta s$')
    plt.ylabel('Norm')
    plt.title('Convergence based on varying ' + r'$\Delta s$' + ' and ' + r'$\Delta t$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, np.round(ds, 3) ))
    dt_values_str = '_'.join(map(str, np.round(dt, 3)))

    # Save the plot to a file -- labels with da values and dt 
    plt.savefig('da_plot/' + folder + '/varied_dt/lw-ex_plot_conv_mu_' + str(c) + '_ds_' + ds_values_str + '_dt_' + dt_values_str + '_order_'+ str(order)  +'.png', dpi=300)  
 
    plt.show()

    return Norm2, L2norm, NormMax, LMaxnorm