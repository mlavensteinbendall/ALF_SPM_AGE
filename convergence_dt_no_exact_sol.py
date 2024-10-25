import numpy as np
import matplotlib.pyplot as plt

def convergence_dt_plt(Smax, ds, dt, order):

    Ntest = len(dt)

    Norm2 = np.zeros([Ntest-1])
    NormMax = np.zeros([Ntest-1])

    L2norm = np.zeros([Ntest-2])
    LMaxnorm = np.zeros([Ntest-2])

    for i in range(0, Ntest-1):

        size = np.arange(0, Smax + ds[i], ds[i])  # Create an array for mesh size
        Nsize = len(size)

        # Load in relevant data for both mesh sizes
        data1 = np.loadtxt(f'da_convergence/num_{i}.txt') 
        data2 = np.loadtxt(f'da_convergence/num_{i+1}.txt')

        # Time of interest to compare
        tinterest = 0.5
        n1 = int(tinterest/dt[i])     # Time step index for data1
        n2 = int(tinterest/dt[i+1])   # Time step index for data2

        # Solve for L2 and L-max norms
        Norm2[i] = np.sqrt(np.mean((data1[n1, :] - data2[n2,  ::2])**2))  # L2 norm
        NormMax[i] = np.max(np.abs(data1[n1, :] - data2[n2,  ::2]))       # L∞ norm

        # Calculate the order of convergence for norms
        if i > 0:
            print(ds[i])
            print(ds[i-1])
            L2norm[i-1] = np.log(Norm2[i] / Norm2[i-1]) / np.log(ds[i] / ds[i-1])
            LMaxnorm[i-1] = np.log(NormMax[i] / NormMax[i-1]) / np.log(ds[i] / ds[i-1])

    # Display the norms and errors
    for i in range(0, Ntest-1):
        print(f'For dt = {round(ds[i], 10)} vs dt = {round(ds[i+1], 10)}')
        print(f'Norm 2 : {round(Norm2[i], 10)}')
        print(f'Norm inf : {round(NormMax[i], 10)}')

        if i > 0:
            print(f'L2 q error: {round(L2norm[i-1], 10)}')  # L2 order of convergence
            print(f'LMax q error: {round(LMaxnorm[i-1], 10)}')  # L∞ order of convergence
        print('')

    # Plot the log-log for the errors
    plt.loglog(dt[:-1], Norm2, label='Norm2')
    plt.loglog(dt[:-1], NormMax, label='NormMax')
    plt.loglog(dt[:-1], dt[:-1]**order, label=f'Order-{order}')

    plt.xlabel('dt')
    plt.ylabel('Norm')
    plt.title('Convergence based on dt')
    plt.legend()
    plt.show()

    return Norm2, L2norm, NormMax, LMaxnorm

# Test if convergence is working
# Mesh options and dt
# da = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
# dt = 0.5 * da  # Time steps based on mesh size

# # Run the function with these parameters
# Smax = 30.0  # Example Smax
# order = 2   # Example order of accuracy

# convergence_dt_plt(Smax, da, dt, order)