import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def convergence_dt_plt(Smax, Tmax, ds, dt, order, c):

    Ntest = len(ds)
    print(Ntest)

    Norm2 = np.zeros([Ntest])
    NormMax = np.zeros([Ntest])

    L2norm = np.zeros([Ntest])
    LMaxnorm = np.zeros([Ntest])


    for i in range(0, Ntest):

        Nage = int(Smax/ds[i])+1
        age = np.zeros([Nage])
        age = np.linspace(0, Smax, Nage) # Create an array of those sizes.

        n = int(Tmax/dt[i]) + 1 # Time-step of comparison.
        Tend = n*dt[i] # Get the associated timepoint value.

        data = np.zeros([int(Tmax/ds[i])+1, Nage])
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
        NormMax[i]  = np.max( np.abs( data[n-1,:] - sol[:] ) )                         # L-Max error.


    # Calculates the L norms -- comparing with the last (Note: ds is increasing)
    for ii in range(0, Ntest - 1):
        L2norm[ii+1]    = np.log( Norm2[ii+1]   / Norm2[ii] )   / np.log( dt[ii+1] / dt[ii] )
        LMaxnorm[ii+1]  = np.log( NormMax[ii+1] / NormMax[ii] ) / np.log( dt[ii+1] / dt[ii] )



    for i in range(0, Ntest):

        print('For ds ='    + str( round( ds[i],10      ) ) + ' and dt ='    + str( round( dt[i],10      ) ) )
        print('Norm 2 error: '   + str( round( Norm2[i], 10  ) ) )
        print('Norm inf error: ' + str( round( NormMax[i], 10) ) )

        if i > 0:
            print('L2 q order: '     + str( round( L2norm[i-1]   , 10    ))) # L2 q estimate.
            print('LMax q order: '   + str( round( LMaxnorm[i-1] , 10    ))) # L-Max q estimate.
            print(' ')

    # plt.figure(figsize=(8, 4))

    # Plot the log-log for the errors.
    plt.loglog(ds, Norm2, label='Norm2')
    plt.loglog(ds, NormMax, label='NormMax')
    plt.loglog(ds, ds**(order), label=f'order-{order }')
    # plt.loglog(ds, ds**1, label=f'order-{1 }')

    # plt.loglog(dt, Norm2, label='Norm2')
    # plt.loglog(dt, NormMax, label='NormMax')
    # plt.loglog(dt, dt**1, label=f'order-{1 }')


    plt.xlabel(r'$\Delta s$')
    plt.ylabel('Norm')
    plt.title('Convergence based on varying ' + r'$\Delta s$' + ' and ' + r'$\Delta t$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, np.round(ds, 3) ))
    dt_values_str = '_'.join(map(str, np.round(dt, 3)))

    # Save the plot to a file -- labels with da values and dt 
    # plt.savefig('ds_plot/varied_dt/plot_conv_mu_' + str(c) + '_ds_' + ds_values_str + '_dt_' + dt_values_str + '_order_'+ str(order)  +'.png', dpi=300)  
 
    plt.show()

    combine = [Norm2, L2norm, NormMax, LMaxnorm]

    return combine