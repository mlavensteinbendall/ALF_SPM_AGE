import numpy as np
import matplotlib.pyplot as plt
from function_trapezoidal_rule import trapezoidal_rule

def reproduction(N, age, da, rep, constant):
    """Calculate the reproduction 
    
    Args:
        N     (array): the number of individuals of each age
    
    Returns:
        births (float): number of births
    """

    if constant:
        reporoduction_rate = np.full(len(age), rep)

    else:
        reporoduction_rate = rep * age

    births = trapezoidal_rule(reporoduction_rate * N, da)

    # print("type : " + str( type(births)))

    return births


import numpy as np

def k_ind(a, u, age_max):
    """
    Function to calculate k, the reproduction rate.
    Daphnia Manga can have clutches up to 100 eggs every 3-4 days until death.
    """
    # Parameters
    p = 2.0
    q = 400
    b = 8.0 / 3.0
    
    # Calculate total population using the helper function
    total_pop = total_population(u, age_max)
    
    # Return k value
    k = (b * np.exp(-a / 10.0) * q**p) / (q**p + total_pop**p)
    return k

def total_population(u, age_max):
    """
    Calculate the total population for a given time-step.
    """
    # Sum up the population values for all age classes
    total_pop = sum(u[:age_max])
    return total_pop

def u_t0(u, age_max):
    """
    Compute the integral using the trapezoidal rule to find the boundary condition.
    """
    s = 0  # Initialize sum for reproduction
    
    # Sum from age 8 to age_max (0-based indexing means starting at 7)
    for i in range(7, age_max):
        # The (i+1) for k_ind(i+1, ...) simulates the 'u[i+1]' aspect in the time step
        s += k_ind(i + 1, u, age_max) * u[i]
        
    return s
