"""
Calculate the time at which the reaction 
stops for a given initial concentration of X.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple
import palettable as pl
import cmasher as cmr
from chemicalclock import ChemicalClock
from helper import find_stop

# Do the usual stuff.
init_conc = lambda a: np.array([1, 0, a, 0, 0, 10, 0, 0])

# Define the time points.
tau = np.linspace(0, 40, 1000)

# Define the rate coefficients.
rate_factor = 1
fixed = 0.1
second_par = fixed/rate_factor
rates = np.array([rate_factor, second_par])

solutions = []
for i in np.linspace(0.001, 10, 100):
    # Create the chemical clock.
    clock = ChemicalClock(init_conc(i), tau, rates, p_fast_factor=100, q_fast_factor=100)

    # Solve the model equations.
    sol = clock.solve()

    # Normalize the solutions.
    sol.y = sol.y / np.max(sol.y, axis=1)[:, None]

    # Find the stop time.
    stop_time = find_stop(sol.y[2])

    # Plot the results.
    solutions.append(stop_time)

# Plot the results.
fig, ax = plt.subplots(1, 1)
ax.plot(np.linspace(0.001, 10, 100), solutions, c="#97e364")
ax.set_xlabel(r"Initial concentration of $x$")
ax.set_ylabel(r"Stop time")
plt.title(r"Stop time as a function of initial concentration of $x$")
ax.grid(color="#424656", alpha=0.1)
plt.show()

