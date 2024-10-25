import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

from function_mortality import mortality

def convergence_dt_plt(Smax, Tmax, ds, dt, order, m, b, constant, folder):
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

        mu = np.zeros(Nstep)
        mu = mortality(Smax, step, m, b, constant)
        # print(mu)

        Ntime = int(Tmax/dt[i])     # Time-step of comparison.
        Tend = Ntime * dt[i]        # Get the associated timepoint value.

        # data = np.zeros([int(Tmax/ds[i]), Nstep])       # initialize matrix for numerical solution
        data = np.zeros([Ntime, Nstep])                   # initialize matrix for numerical solution
        sol = np.zeros([Nstep])                           # initialize array for analytical solution


        # Numerical solution -- download relevent data
        data = np.loadtxt('da_convergence/num_' + str(i) + '.txt') 

        # Analyticial solution -- changes for what ds is
        # sol = np.exp(-(step - ( Tend + 5))**2) * np.exp(-mu * Tend)     # with advection -- constant
        # sol = np.exp(-(step - ( Tend + 5))**2) * np.exp( - m *  step  * Tend + m * Tend**2 * 0.5)     # with advection -- non-constant
        # sol = np.exp(-(step - 5)**2) * np.exp( - mu * Tend)             # without advection

        sol = np.exp(-(step - ( Tend + 5))**2) * np.exp(- m * step * Tend + 0.5 * m * (Tend)**2)     # with advection -- NON CONSTANT - slope

        # sol = np.exp(-(step - ( Tend + 5))**2) * np.exp(- (Tend + 200 * np.log(step**2 + 400) - 200 * np.log((step - Tend)**2 +400))) # hill function
        # sol = np.exp(-(step - ( Tend + 5))**2) * np.exp(- (20 * np.log(step**2 + 400) - 20 * np.log((step - Tend)**2 +400)))  # hill function


    

        # Calculate the norm 
        Norm2[i] = np.sqrt(np.mean((data[-1, :] - sol[:])**2))      # L2 norm
        NormMax[i]  = np.max( np.abs( data[-1,:] - sol[:] ) )       # Lmax norm


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

    plt.loglog(dt, Norm2, label='Norm2')
    plt.loglog(dt, NormMax, label='NormMax')
    plt.loglog(dt, dt**(order), label=f'order-{order }')

    plt.xlabel(r'$\Delta t$')
    plt.ylabel('Norm')
    plt.title('Convergence based on varying ' + r'$\Delta s$' + ' and ' + r'$\Delta t$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, np.round(ds, 3) ))
    dt_values_str = '_'.join(map(str, np.round(dt, 3)))

    # Save the plot to a file -- labels with da values and dt 
    # plt.savefig('da_plot/' + folder + '/varied_dt/lw-ex_plot_conv_mu_' + str(c) + '_ds_' + ds_values_str + '_dt_' + dt_values_str + '_order_'+ str(order)  +'.png', dpi=300)  
 
    plt.show()

    return Norm2, L2norm, NormMax, LMaxnorm

# Test function
# Mesh options and dt
# da = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
# dt = 0.5 * da  # Time steps based on mesh size

# # Run the function with these parameters
# Smax = 30.0  # Example Smax
# Tmax = 5
# order = 2   # Example order of accuracy
# m = 0.5
# b = 0
# constant = False

# convergence_dt_plt(Smax, Tmax, da, dt, order, m, b, constant, 'folder')
