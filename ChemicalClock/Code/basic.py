"""
Try to get a basic plot of oscillations in
the concentrations of Iodine and Iodide.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import palettable as pl
import cmasher as cmr
from chemicalclock import ChemicalClock

# Define the initial concentrations.
init_conc = np.array([1, 1, 10, 0, 0, 10, 0, 0])

# Define the time points.
tau = np.linspace(0, 100, 1000)

# Define the rate coefficients.
rates = np.array([0.1, 1])
# p_fast_factor = 1
# q_fast_factor = 1

# Create the chemical clock.
clock = ChemicalClock(init_conc, tau, rates, p_fast_factor=2, q_fast_factor=10)

# Solve the model equations.
sol = clock.solve()

# Normalize to relative concentrations.
# sol.y = sol.y / np.sum(init_conc)

# Plot the results.
fig, ax = plt.subplots(1, 1)
ax.plot(sol.t, sol.y[0], label='[M]')
ax.plot(sol.t, sol.y[1], label='[N]')
# ax.plot(sol.t, sol.y[2], label='[U]')
# ax.plot(sol.t, sol.y[3], label='[V]')
# ax.plot(sol.t, sol.y[4], label='[W]')
# ax.plot(sol.t, sol.y[5], label='[X]')
# ax.plot(sol.t, sol.y[6], label='[Y]')
# ax.plot(sol.t, sol.y[7], label='[Z]')
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Relative Concentration')
ax.set_title('Full solution for a Iodine Chemical Clock')
plt.show()
