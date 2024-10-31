import numpy as np
import matplotlib.pyplot as plt

def convergence_da_plt(Smax, da, dt, order, folder):

    print('Calculate convergence varying da and fixing dt (no analytic solution)')

    Ntest    = len(da)

    Norm2    = np.zeros([Ntest-1])
    NormMax  = np.zeros([Ntest-1])

    L2norm   = np.zeros([Ntest-2])
    LMaxnorm = np.zeros([Ntest-2])

    for i in range(0, Ntest-1):

        # Load in relevant data for both mesh sizes
        data1 = np.loadtxt(f'da_convergence/num_{i}.txt') 
        data2 = np.loadtxt(f'da_convergence/num_{i+1}.txt')

        # Time of interest to compare
        tinterest = 2.5
        n = int(tinterest/dt)     # Time step index for data1
        # n2 = int(tinterest/da[i+1])   # Time step index for data2

        # Solve for L2 and L-max norms
        Norm2[i] = np.sqrt(np.mean((data1[n, :] - data2[n,  ::2])**2))  # L2 norm
        NormMax[i] = np.max(np.abs(data1[n, :]  - data2[n,  ::2]))       # L∞ norm

        # Calculate the order of convergence for norms
        if i > 0:
            print(da[i])
            print(da[i-1])
            L2norm[i-1]   = np.log(Norm2[i-1]   / Norm2[i])   / np.log(da[i-1] / da[i])
            LMaxnorm[i-1] = np.log(NormMax[i-1] / NormMax[i]) / np.log(da[i-1] / da[i])

    # Display the norms and errors
    for i in range(0, Ntest-1):
        print(f'For da =   {round(da[i], 10)} vs da = {round(da[i+1], 10)}')
        print(f'Norm 2 :   {round(Norm2[i], 10)}')
        print(f'Norm inf : {round(NormMax[i], 10)}')

        if i > 0:
            print(f'L2 q error:   {round(L2norm[i-1], 10)}')  # L2 order of convergence
            print(f'LMax q error: {round(LMaxnorm[i-1], 10)}')  # L∞ order of convergence
        print('')

    # Plot the log-log for the errors
    plt.loglog(da[:-1], Norm2,          label='Norm2')
    plt.loglog(da[:-1], NormMax,        label='NormMax')
    plt.loglog(da[:-1], da[:-1]**order, label=f'Order-{order}')

    plt.xlabel(r'$\Delta a$')
    plt.ylabel('Norm')
    plt.title('Convergence based on varying ' + r'$\Delta a$' + ' and fixing' + r'$\Delta t$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, np.round(da, 3) ))
    dt_values_str = dt


    # Save the plot to a file -- labels with da values and dt 
    # plt.savefig('da_plot/' + folder + '/varied_dt/lw-ex_plot_conv_mu_' + str(0) + '_da_' + ds_values_str + '_dt_' + dt_values_str + '_order_'+ str(order)  +'.png', dpi=300)  
    plt.show()

    return Norm2, L2norm, NormMax, LMaxnorm


# da = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
# dt = 0.0001

# # # Run the function with these parameters
# Smax = 30.0  # Example Smax
# order = 2   # Example order of accuracy

# convergence_da_plt(Smax, da, dt, order, "test")
