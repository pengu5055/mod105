"""
Basic test and run of the BinaryEquilibrium class. So solving
binary reactions in the approximation of equilibrium.
"""
import numpy 
import matplotlib.pyplot as plt
import palettable as pl
from binaryreaction import BinaryEquilibrium

# Define the initial concentrations
k = 1000
s = 100

# Define the initial concentrations
init_conc = [1, 0, 0, 0]

# Define the time at which to evaluate the solution
tau = numpy.linspace(0, 10, 1000)

# Initialize the BinaryReaction object
br = BinaryEquilibrium(init_conc, tau, k, s)

# Solve the system
sol = br.solve()

# Plot the results
colors = pl.colorbrewer.sequential.Purples_7.hex_colors

plt.plot(tau, sol[0], c=colors[2], label='a')
plt.plot(tau, sol[1], c=colors[3], label='a*')
plt.plot(tau, sol[2], c=colors[4], label='b')
plt.plot(tau, sol[3], c=colors[5], label='c')

plt.legend()
plt.show()
