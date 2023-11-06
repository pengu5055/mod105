"""
The idea here is to plot a matrix depicting the rate of change
of the concentration of a solved for different values of k and s.
The rate of change should be obtained by fitting an exponential
function to the solution and then taking the decay constant as
the rate of change.
"""
import numpy as np
import matplotlib.pyplot as plt
import palettable as pl
from binaryreaction import BinarySingular
from helper import *

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

colors = pl.cartocolors.sequential.agSunset_7.hex_colors
colors2 = pl.cartocolors.sequential.TealGrn_7.hex_colors

plt.plot(tau, sol[0], c=colors[0], label='a')
fit, par, err = find_decay_parameter(tau, sol[0], return_curve_fit=True)

plt.plot(tau, fit, c=colors2[0], label='fit')

# Create textbox with decay parameter

textstr = '\n'.join((
    r'$\tau=%.2f$' % (par[1], ),
    r'$\sigma=%.2f$' % (err[1, 1], )))
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

plt.text(0.84, 0.14, textstr, transform=plt.gca().transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.legend()
plt.show()
