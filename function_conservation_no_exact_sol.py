import numpy as np
import matplotlib.pyplot as plt

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

def conservation_plt(Smax, da, dt, order, folder):

    print('Calculate conservation (no analytic solution)')

    Ntest    = len(da)

    Norm1    = np.zeros([Ntest-1])

    L1norm   = np.zeros([Ntest-2])

    totalPop_1 = np.zeros([5]) 
    totalPop_2 = np.zeros([5]) 

    for i in range(0, Ntest-1):

        # Load in relevant data for both mesh sizes
        data1 = np.loadtxt(f'da_convergence/num_{i}.txt') 
        data2 = np.loadtxt(f'da_convergence/num_{i+1}.txt')

        # Time of interest to compare
        # tinterest = 2.5
        # n1 = int(tinterest/da[i])     # Time step index for data1
        # n2 = int(tinterest/da[i+1])   # Time step index for data2

        totalPop_1[i] = trapezoidal_rule( data1[-1,:], da[i])           # using data at Tmax
        totalPop_2[i] = trapezoidal_rule( data2[-1,:], da[i+1])         # using data at Tmax
        print('Numerical total pop using ' + r'$\Delta a$' + ' = ' + str(da[i]) +'  = '      + str(totalPop_1[i]))
        print('Numerical total pop using ' + r'$\Delta a$' + ' = '  + str(da[i+1]) +'  = '    + str(totalPop_2[i]))

        # Solve for L2 and L-max norms
        Norm1[i] = np.abs( totalPop_1[i] - totalPop_2[i])  # L2 norm


        # Calculate the order of convergence for norms
        if i > 0:
            L1norm[i-1]   = np.log(Norm1[i-1]   / Norm1[i])   / np.log(da[i-1] / da[i])

    # Display the norms and errors
    for i in range(0, Ntest-1):
        print(f'For da =   {round(da[i], 10)} vs da = {round(da[i+1], 10)}')
        print(f'Norm 1 :   {Norm1[i]}')


        if i > 0:
            print(f'L1 q error:   {L1norm[i-1]}')  # L2 order of convergence
        print('')

    # Plot the log-log for the errors
    plt.loglog(da[:-1], Norm1,          label='Norm1')
    plt.loglog(da[:-1], da[:-1]**order, label=f'Order-{order}')

    plt.xlabel(r'$\Delta a$')
    plt.ylabel('Norm')
    plt.title('Convergence based on varying ' + r'$\Delta a$' + ' and fixing' + r'$\Delta t$')
    plt.legend()

    # Convert ds array values to a string
    ds_values_str = '_'.join(map(str, np.round(da, 3) ))

    # if isinstance(dt, np.ndarray):
    if isinstance(dt, np.ndarray):
        plt.savefig('da_plot/'+ folder +'/varied_dt/lw-ex_plot_totPop_mu__ds_' + ds_values_str + '.png', dpi=300)  

    else:
        plt.savefig('da_plot/'+ folder +'/fixed_dt/lw-ex_plot_totPop_mu__ds_' + ds_values_str + '.png', dpi=300)  
 
    plt.show()

    return Norm1, L1norm


# da = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
# dt = 0.0001

# # # Run the function with these parameters
# Smax = 30.0  # Example Smax
# order = 2   # Example order of accuracy

# conservation_plt(Smax, da, dt, order, "test")
