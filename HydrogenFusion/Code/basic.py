"""
Basic evaluation of the Hydrogen Fusion process.
"""
import numpy as np
import matplotlib.pyplot as plt
import palettable as pl
from hydrogenfusion import HydrogenFusion

# Define the initial concentrations
init_conc = np.array([1, 1, 0, 0, 0])

# Define the rate coefficients
rates = np.array([1, 1, 1, 1, 1])

# Define the time points to evaluate the rate equations
t_eval = np.linspace(0, 10, 1000)

# Initialize the HydrogenFusion class
hf = HydrogenFusion(init_conc, t_eval, rates)

# Solve the rate equations
sol = hf.solve()

# Plot the results
fig, ax = plt.subplots()

ax.plot(sol.t, sol.y[0], label='u')
ax.plot(sol.t, sol.y[1], label='v')
ax.plot(sol.t, sol.y[2], label='x')
ax.plot(sol.t, sol.y[3], label='y')
ax.plot(sol.t, sol.y[4], label='z')

ax.set_xlabel('Time (s)')
ax.set_ylabel('Concentration (M)')

ax.legend()
plt.show()

