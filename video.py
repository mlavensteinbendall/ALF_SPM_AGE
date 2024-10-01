import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the relevant data
population_data = np.loadtxt('da_convergence/num_4.txt') 

age_max = 15
da = 0.00625

# Number of time steps is the number of outer elements
time_steps = len(population_data)

# Number of age groups is the number of inner elements
age_groups = len(population_data[0])

# Set up the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(0, age_max)  # Rescale the x-axis to be from 0 to 15
ax.set_ylim(0, np.max(population_data))  # Population numbers on the y-axis
line, = ax.plot([], [], lw=2)

# Create an array of age labels from 0 to age_max
age_labels = np.linspace(0, age_max, age_groups)

# Initialize the plot with empty data
def init():
    line.set_data([], [])
    return line,

# Update function that will be called at each time step
def update(frame):
    # X-axis: Corresponding age labels, Y-axis: Population numbers at the current time step
    y = population_data[frame]
    line.set_data(age_labels, y)  # Use the age labels for the x-values
    ax.set_title(f"Distribution of Population at Time : {round(10 * frame * 0.00625 * 0.5, 0)}", fontsize=16)  # Update the title
    return line,

# Create the animation
ani = FuncAnimation(fig, update, frames=np.arange(0, time_steps), init_func=init, blit=False, interval=50)

# Set axis labels
plt.xlabel('Age', fontsize=14)
plt.ylabel('Distribution of Population', fontsize=14)

# Set tick label sizes without changing tick positions
ax.tick_params(axis='x', labelsize=12)  # Set x-axis tick label size
ax.tick_params(axis='y', labelsize=12)  # Set y-axis tick label size


# Save the animation
ani.save('population_dynamics_mu_0.gif', writer='pillow', fps=30)  # Save as GIF
# ani.save('population_dynamics.mp4', writer='ffmpeg', fps=30)  # Uncomment to save as MP4


# Show the animation
plt.show()