"""
A test file to try out in development code.
"""
import numpy as np
import matplotlib.pyplot as plt
import palettable as pl
from binaryreaction import BinarySingular

# Define the initial concentrations
init_conc = [1, 0, 0, 0]

# Define the time at which to evaluate the solution
tau = np.linspace(0, 10, 1000)

# Define the rate constants
k = 1000
s = 0.1

# Initialize the BinaryReaction object
br = BinarySingular(init_conc, tau, k, s)

# Solve the system
sol = br.solve()

# Plot the results

colors = pl.cartocolors.sequential.TealGrn_7.hex_colors

plt.plot(tau, sol[0], c=colors[0], label='a')
plt.plot(tau, sol[1], ls=":", c=colors[1], label='a*')
plt.plot(tau, sol[2], ls="-", c=colors[2], label='b')
plt.plot(tau, sol[3], ls="-", c=colors[3], label='c')

plt.legend()
plt.show()
