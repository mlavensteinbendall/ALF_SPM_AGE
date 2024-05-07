
import numpy as np # Numpy for numpy
import matplotlib.pyplot as plt

def LW_SPM(ds,dt,ntag,filename):

    Smax = 15 # Maximum size to calculate
    Tmax = 1 # End time
    Nsizes = int(Smax/ds)+1 # Total number of size-steps
    Ntimes = int(Tmax/dt)+1 # Total number of time-steps
    # sizes = np.linspace(ds,Smax,num=Nsizes) # Grid of size points
    sizes = np.linspace(0,Smax,num=Nsizes) # Grid of size points
    times = np.linspace(0,Tmax,num=Ntimes) # Grid of times points

    # Initial condition
    N = np.zeros([Nsizes]) # Store for the current timepoint

    # Fitness functions
    g = np.ones([Nsizes]) # growth rate
    g[:] = np.exp(-(sizes[:])) # To do: test for exp(- size/6)
    r = np.zeros([Nsizes])  # reproduction
    r[int(Nsizes/4):-1] = 0
    mu = np.zeros([Nsizes]) # mortality
    # mu = np.ones([Nsizes]) # mortality
    # mu = sizes #  To do: Look at smaller mu x


    #Plot the growth rate and mortality rate
    # plt.plot(sizes,g)
    # plt.plot(sizes,mu)
    # plt.show()

 

    # Difference matrices
    D1 = np.zeros([Nsizes,Nsizes]) # Store for 1st finite difference matrix
    D2 = np.zeros([Nsizes,Nsizes]) # Store for 2nd finite difference matrix

    # End-points
    # D1[0,0]  = -3/(2*ds); D1[0,1]    = 2/(ds);  D1[0,2] = -1/(2*ds)
    # D1[-1,-1] = 3/(2*ds); D1[-1,-2] = -2/(ds); D1[-1,-3] = 1/(2*ds)
    # D2[0,0]  = 2/(ds**2); D2[0,1] = -5/(ds**2);  D2[0,2] = 4/(ds**2); D2[0,3] = -1/(ds**2)
    # D2[-1,-1] = 2/(ds**2); D2[-1,-2] = -5/(ds**2); D2[-1,-3] = 4/(ds**2); D2[-1,-4] = -1/(ds**2)

    # Mid points
    for ii in range(1,Nsizes-2): # Loop through the diagonals
        D1[ii,ii-1] = -0.5/ds;   D1[ii,ii+1] = 0.5/ds  # 1st derivative centeral finite-difference
        D2[ii,ii-1] = 1/(ds**2); D2[ii,ii] = -2/(ds**2); D2[ii,ii+1] = 1/(ds**2) # 2nd derivative centeral difference 
        # D1[ii,ii-2] = (1/12)/ds; D1[ii,ii-1] = (-2/3)/ds; D1[ii,ii+1] = (2/3)/ds; D1[ii,ii+2] = (-1/12)/ds
        # D2[ii,ii-2] = (-1/12)/(ds**2); D2[ii,ii-1] = (4/3)/(ds**2); D2[ii,ii] = (-5/2)/(ds**2); D2[ii,ii+1] = (4/3)/(ds**2); D2[ii,ii+2] = (-1/12)/(ds**2)


    a1 = np.zeros([Nsizes]) # Store for N terms coef
    a2 = np.zeros([Nsizes]) # Store for dNds terms coef
    a3 = np.zeros([Nsizes]) # Store for d2Nds2 terms coef

    a1[:] = 1
    a2[:] = -1*dt 
    a3[:] = (dt**2)/2

    # Initial condition
    # N[:] = np.exp(-(sizes-10)**2) # Initial condition setter  -- Gaussian centered at 10
    # N[:] = np.exp(-((sizes-2)/0.1)**2) # change more to the right
    # N[:] = np.exp(-((sizes-0.4)/0.1)**2)
    N[:] = np.exp(-((sizes-5)/1)**2) #need to be wide enought to fix the oscilations 
    # N[:] = np.exp(-((sizes-2)/0.5)**2)

    # CFL dt < exp(-0.4) ds
    # divide by 2 10 something smaller t han 1

    # plt.plot(sizes, N)
    # plt.xlabel('Size')
    # plt.ylabel('Population')
    # plt.title('Population Based on Size at Initial Time')
    # plt.show()  

    with open(filename + 'upwind_num_' + str(ntag) + '.txt', 'w') as file: # Initialise an outputter file (safe)
        for t,T in enumerate(times): # Loop on times

            for n in N: # Output the current time solution
                file.write(str(n))
                file.write(" ")
            file.write("\n")

            N[0] = 0

            # Step 1 - half step time, mortality
            N[:] = N[:]*np.exp(-mu[:]*dt/2)

            # Step 2 - half step time, growth
            N[:] = a1[:]*N[:] + a2[:]*(D1.dot(N)) + a3[:]*(D2.dot(N))

            # Step 3 - half step time, mortality
            N[:] = N[:]*np.exp(-mu[:]*dt/2)

            
        for n in N: # Output the final time solution
            file.write(str(n))
            file.write(" ")
        file.write("\n")
            
    # # Print the solution over time
    # plt.plot(N)
    # plt.title('N_final')
    # plt.show()

    # size_values = np.arange(len(N)) * ds

    # plt.plot(size_values, N)
    # plt.xlabel('Size')
    # plt.ylabel('Population')
    # plt.title('Population Based on Size at Final Time for ds = ' + str(ds))
    # plt.show()

    
    # plt.plot(sizes, N)
    # # plt.title('N0')
    # plt.show()      
        

   # Check to see if there are negative here
    # result = a1[:] - 2*a3[:]
    # result = a3[:] - a2[:]

    # # Check if there are any negative values
    # negative_values_exist = any(value < 0 for value in result)

    # if negative_values_exist:
    #     print("There are negative values in the expression (a1[:] - 2*a3[:])")
    # else:
    #     print("There are no negative values in the expression (a1[:] - 2*a3[:])")