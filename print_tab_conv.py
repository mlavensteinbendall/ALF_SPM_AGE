from tabulate import tabulate
import numpy as np

def tabulate_conv(dt, ds, Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm):
    # Define the headers
    headers = ["$\Delta{t}$","$\Delta{s}$", "$||N-N_{ex}||_2$", "$q_2$", "$||N-N_{ex}||_\infty$", "$q_\infty$", "$||N-N_{ex}||_1$", "$q_1$"]

    time = np.zeros([len(ds)])
    time[:] = dt

    # Combine the data into a list of tuples
    latex = list(zip(time, ds,  Norm2, L2norm, NormMax, LMaxnorm, Norm1, L1norm))

    # Generate the LaTeX table
    latex_table = tabulate(latex, headers=headers, tablefmt="latex_raw")

    # Print or save the LaTeX table
    print(latex_table)
