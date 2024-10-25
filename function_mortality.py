import numpy as np
import matplotlib.pyplot as plt

def mortality(age_max, age, m, b, constant):
    """Calculates the numerical solution using strang splitting, lax-wendroff, and runge-kutta method. 
    
    Args:
        age     (array): A list of all the ages
        m       (int):   A constant for the slope
        b       (int):   y-intercept
        constant(bool):  Check whether a constant mortality is wanted or not
        
    Returns:
        mu       (array): Represents mortality rate at each age
    """
    if constant:
        # Apply constant mortality rate
        mu = np.full(len(age), m)  # Fill array with constant value of m (assuming m in [0,1])

    else:
        # Apply mortality based on the linear equation: y = m * (age / age_max) + b
        # mu = m *(1 - np.cos(age/age_max * np.pi))
        mu = m * age + b

        # mu = 1 + age * (20**2 / (20**2 + age**2))
        # mu = age/10 * (20**2 / (20**2 + age**2))

        # mu = 0.1 + 0.9 * np.exp(- 1000000000 * np.exp(- 1 * age))

        # mu = m * (age/age_max) + b
        # mu = m * (age)/ (age+10000) + b
        # mu = m * (age)/ (age+1000) + b
        # mu =  np.exp( - age )

        # mu = np.log(age+1)


    # # clear figures 
    # plt.clf()

    # # Plot.
    # plt.plot(age, mu)

    # plt.xlabel('Age')
    # plt.ylabel('Mortality Rate')
    # plt.title('Mortality Rate based on Age')
    # plt.legend()
    # plt.show()

    return mu


