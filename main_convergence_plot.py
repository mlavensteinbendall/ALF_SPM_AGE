import numpy as np
# from convergence_da import convergence_da_plt
from convergence_dt import convergence_dt_plt
from function_conservation import conservation_plt
from print_tab_conv import tabulate_conv
import matplotlib.pyplot            as plt 

## INITIAL CONDITIONS
Amax = 30       # max age
# Tmax = 5     # max time
Tmax = 20
order = 2       # order of method
Ntest = 5       # number of cases


# Mortality set up
m = 0 #1/30            # constant for mux
b = 0              # y-intercept

# Reproduction set up
k = 1

# need to chose da and dt so that the last value in the array are Amax and Tmax
da = np.zeros([Ntest]) # order smallest to largest

# # vary da and dt cases:
da[0] = 0.1
da[1] = 0.05
da[2] = 0.025
da[3] = 0.0125
da[4] = 0.00625

dt = np.zeros([Ntest]) # order smallest to largest

# dt[0] = 0.5 * 0.1
# dt[1] = 0.5 * 0.05
# dt[2] = 0.5 * 0.025
# dt[3] = 0.5 * 0.0125
# dt[4] = 0.5 * 0.00625

dt = 0.5 * da
# dt = 0.01
# dt = 0.001
# dt = 0.0001

time = np.arange(0,Tmax + dt[0], dt[0])
age = np.arange(0, Amax + da[0], da[0])

testing_folder = "mortality"
function_folder = "constant_mortality"

if isinstance(dt, np.ndarray):
    convergence_folder = 'varied_dt'

else:
    convergence_folder = 'fixed_dt'


# folder = 'convergence/' + testing_folder + '/' + function_folder + '/' + convergence_folder

# Norm2, L2norm, NormMax, LMaxnorm = convergence_dt_plt(Amax, Tmax, da, dt, 2, m, b, True, folder)

# Norm1, L1norm = conservation_plt(Ntest, da, m, Amax, Tmax, dt, 2, folder, True, False)

# tabulate_conv(dt, da, Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm, folder)

print('Plot boundary condition')

folder = '/Users/mlavensteinbendall/Documents/ALF_SPM_AGE/convergence/reproduction/no_mortality/constant_reproduction/rep_1-original/varied_dt/solutions/'

data = np.loadtxt(f'{folder}/num_0.txt') 

print(data[55,0])
print(data[55,1])
print(data[55,2])

print(data[100,0])
print(data[100,1])
print(data[100,2])
print(data[101,0])
print(data[101,1])
print(data[101,2])

# 55 

# plt.plot(age, data[55,:])
plt.plot(time, data[:,0], label = 'age = 0')
plt.plot(time, data[:,1], label = 'age = da')
plt.plot(time, data[:,2], label = 'age = 2 * da')
plt.plot(time, data[:,3], label = 'age = 3 * da')
plt.plot(time, data[:,4], label = 'age = 4 * da')
plt.plot(time, data[:,5], label = 'age = 5 * da')
plt.plot(time, data[:,40], label = 'age = 40 * da')
plt.plot(time, data[:,60], label = 'age = 60 * da')
plt.plot(time, data[:,100], label = 'age = 100 * da')
plt.plot(time, data[:,250], label = 'age = 250 * da')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Population for different ages over time')
plt.legend()
plt.show()

    # if isinstance(dt, np.ndarray):
    #     plt.savefig(folder + '/plots/boundary_condition_for_da_' + str(da) + '_dt_' + str(dt[index]) + '.png', dpi=300)

    # else:
    #     plt.savefig(folder + '/plots/dt_' + str(dt) + '/boundary_condition_for_da_' + str(da) + '_dt_' + str(dt) + '.png', dpi=300)

    # plt.close()