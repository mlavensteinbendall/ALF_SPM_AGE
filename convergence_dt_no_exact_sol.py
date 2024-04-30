import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def convergence_dt_plt(Smax, Tmax, ds, dt):

    Ntest = len(dt)

    Norm2 = np.zeros([Ntest-1])
    NormMax = np.zeros([Ntest-1])

    L2norm = np.zeros([Ntest-2])
    LMaxnorm = np.zeros([Ntest-2])


    for i in range(0, Ntest-1):

        Nsize   = int(Smax/ds[i+1])+1
        size    = np.zeros([Nsize])
        size    = np.linspace(0,Smax,Nsize) # Create an array of those sizes.

        data1   = np.loadtxt('ds_convergence/upwind_num_' + str(i) + '.txt') # Load in relevant data.
        data2   = np.loadtxt('ds_convergence/upwind_num_' + str(i+1) + '.txt') # Load in relevant data.

        # Calculate the equivilent 
        tinterest = 0.5 # Define the time where the calculations are compared.
        n1      = int(tinterest/dt[i]) # Get fine-data time-step.
        n2      = int(tinterest/dt[i+1]) # Get corse-data time-step.

        # Solve for L-2 and L-max
        Norm2[i]    = ( ( 1 / len(size) ) * np.sum( ( data1[n1,0::2] - data2[n2,:] ) **2 ) ) **0.5 # L2 error.
        NormMax[i]  = np.max( np.abs( data1[n1,0::2] - data2[n2,:] ) )                             # L-Max error.

        if i > 0 :
            # Loop through the remaining datasets.
            for ii in range(0,Ntest-2):
                # L2norm[ii-1]    = np.log( Norm2[ii-1]   / Norm2[ii] )   / np.log( ds[ii-1] / ds[ii] )
                # LMaxnorm[ii-1]  = np.log( NormMax[ii-1] / NormMax[ii] ) / np.log( ds[ii-1] / ds[ii] )

                L2norm[ii]    = np.log( Norm2[ii+1]   / Norm2[ii] )   / np.log( ds[ii+1] / ds[ii] )
                LMaxnorm[ii]  = np.log( NormMax[ii+1] / NormMax[ii] ) / np.log( ds[ii+1] / ds[ii] )


    for i in range(0, Ntest-1):

        print('For dt ='    + str( round( dt[i],10      ) ) + ' v.s. ' 'For dt ='    + str( round( dt[i+1],10      ) ))
        print('Norm 2 : '   + str( round( Norm2[i], 10  ) ) )
        print('Norm inf : ' + str( round( NormMax[i], 10) ) )

        # L2 Norm and # LMax Norm
        print('L2 q error: '     + str( round( L2norm[ii-1]   , 10    ))) # L2 q estimate.
        print('LMax q error: '   + str( round( LMaxnorm[ii-1] , 10    ))) # L-Max q estimate.
        print(' ')

    # Plot the log-log for the errors.
    plt.loglog(Norm2, label='Norm2')
    plt.loglog(NormMax, label='NormMax')
    plt.loglog(dt**2, label='order-2')

    plt.xlabel('dt')
    plt.ylabel('Norm')
    plt.title('Convergence based on dt')
    plt.legend()
    plt.show()