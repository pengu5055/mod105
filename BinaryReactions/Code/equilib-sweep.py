"""
Parameter sweep for BinaryEquilibrium.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmasher as cmr
import palettable as pl
from binaryreaction import BinaryEquilibrium

# Define the initial concentrations
init_conc = [1, 0, 0, 0]

# Define the time at which to evaluate the solution
tau = np.linspace(0, 10, 1000)

# Define the rate constant sweep ranges
k = 1000
s1 = np.linspace(1.1, 10, 10)

# Initialize the BinaryReaction object
solutions = []
solutions2 = []
for par in s1:
    br = BinaryEquilibrium(init_conc, tau, k, par)

    sol = br.solve()

    solutions.append(sol)


# Plot the results
fig, ax = plt.subplots()
cm = pl.colorbrewer.sequential.Purples_7.mpl_colormap
colors = cmr.take_cmap_colors(cm, len(s1))

for i, sol in enumerate(solutions):
    ax.plot(tau, sol[0], c=colors[i], label=f's = {s1[i]:.2f}')

ax.set_xlabel('Dimensionless Time')
ax.set_ylabel('Relative Concentration')
ax.set_title('Equilibrium Approximation')
plt.show()
