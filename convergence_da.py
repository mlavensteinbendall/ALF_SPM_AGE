import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def convergence_da_plt(age_max, time_max, da, dt, order, c):

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

        Nage = int(age_max/da[i])+1
        age = np.zeros([Nage])
        age = np.linspace(0, age_max, Nage) # Create an array of those sizes.

        data = np.zeros([int(time_max/da[i])+1, Nage])
        sol = np.zeros([Nage])

        data = np.loadtxt('ds_convergence/upwind_num_' + str(i) + '.txt') # Load in relevant data.

        # Analyticial solution -- changes for what ds is
        sol = np.exp(-(age - ( Tend + 5))**2) * np.exp(-c * Tend)          # mu(s) = 0
        
        # # plt data vs sol
        # plt.plot(age,data[-1,:]) 
        # plt.plot(age,sol) # looks right
        # plt.show()

        # Solve for L-2 and L-max
        Norm2[i]    = ( ( 1 / Nage ) * np.sum( np.abs( data[n-1,:] - sol[:] ) **2 ) ) **0.5  # L2 error.
        NormMax[i]  = np.max( np.abs( data[n-1,:] - sol[:] ) )                               # L-Max error.


    # Calculates the L norms -- comparing with the last (Note: ds is increasing)
    for ii in range(0, Ntest - 1):
        L2norm[ii+1]    = np.log( Norm2[ii+1]   / Norm2[ii] )   / np.log( da[ii+1] / da[ii] )
        LMaxnorm[ii+1]  = np.log( NormMax[ii+1] / NormMax[ii] ) / np.log( da[ii+1] / da[ii] )



    for i in range(0, Ntest):

        print('For ds ='            + str( round( da[i],        10  ) ) )
        print('Norm 2 error: '      + str( round( Norm2[i],     10  ) ) )
        print('Norm inf error: '    + str( round( NormMax[i],   10  ) ) )

        if i > 0:
            print('L2 q order: '    + str( round( L2norm[i-1]   , 10    ) ) ) # L2 q estimate.
            print('LMax q order: '  + str( round( LMaxnorm[i-1] , 10    ) ) ) # L-Max q estimate.
            print(' ')

    # plt.figure(figsize=(8, 4))

    # Plot the log-log for the errors.
    plt.loglog(da, Norm2, label='Norm2')
    plt.loglog(da, NormMax, label='NormMax')
    plt.loglog(da, da**(order), label=f'order-{order }')
    # plt.loglog(da, da**1, label=f'order-{1 }')


    plt.xlabel(r'$\Delta a$')
    plt.ylabel('Norm')
    plt.title('Convergence based on ' + r'$\Delta a$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, da))

    # Save the plot to a file -- labels with da values and dt 
    # plt.savefig('ds_plot/fixed_dt/plot_conv_mu_'+ str(c) + '_ds_'+ ds_values_str + f'_dt_{dt }'+ '_order_' + str(order) +'.png', dpi=300)  

    plt.show()

    combine = [Norm2, L2norm, NormMax, LMaxnorm]

    return combine

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