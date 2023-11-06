"""
Basic test of the BinaryReactions class.
"""
import numpy as np
import matplotlib.pyplot as plt
import palettable as pl
from binaryreaction import BinaryReaction


# Define the initial concentrations.
init_conc = np.array([1, 0, 0, 0])

# Define the time points.
tau = np.linspace(0, 100, 1000)

# Define the constants.
k = 1
s = 0.1

# Initialize the BinaryReaction object.
br = BinaryReaction(init_conc, tau, k, s)

# Solve the system.
sol = br.solve()

# Plot the results.
colors = pl.colorbrewer.sequential.PuRd_7.hex_colors

plt.plot(tau, sol.y[0], c=colors[2], label='a')
plt.plot(tau, sol.y[1], c=colors[3], label='a*')
plt.plot(tau, sol.y[2], c=colors[4], label='b')
plt.plot(tau, sol.y[3], c=colors[5], label='c')
plt.legend()

plt.ylabel('Relative Concentration')
plt.xlabel('Dimensionless Time')
plt.show()
