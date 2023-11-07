"""
Experimentation with parameters to see if any interesting 
results can be found.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmasher as cmr
import palettable as pl
from hydrogenfusion import HydrogenFusion
from helper import *

# Define the initial concentrations
init_conc = np.array([1, 1, 0, 10, 0])

# Define the rate coefficients
rates = np.array([2, 0.5, 1, 1, 0.1])

# Define the time points to evaluate the rate equations
t_eval = np.linspace(0, 1000, 10000)

# Initialize the HydrogenFusion class
hf = HydrogenFusion(init_conc, t_eval, rates)

# Solve the rate equations
sol = hf.solve()

# Normalize the solutions to the sum of the initial concentrations
# thus giving relative concentrations instead of absolute ones.
sol.y = sol.y / np.sum(init_conc)


# Plot the results
fig, ax = plt.subplots()

cm = pl.colorbrewer.sequential.BuPu_7.mpl_colormap
colors = cmr.take_cmap_colors(cm, 7, cmap_range=(0.2, 1), return_fmt='hex')

# --- Plot 1 ---
ax.plot(sol.t, sol.y[0], label='$H_2$', color=colors[1], linestyle='--')
ax.plot(sol.t, sol.y[1], label='$Br_2$', color=colors[5], linestyle='--')
ax.plot(sol.t, sol.y[2], label='$HBr$', color=colors[3])
ax.plot(sol.t, sol.y[3], label='$H$', color=colors[2], linestyle=':')
ax.plot(sol.t, sol.y[4], label='$Br$', color=colors[4], linestyle=':')

# NOTE: An idea is to make the curves dimensionless by normalizing to the sum
# of the initial concentrations (thus giving relative concentrations).

ax.set_xlabel('Time (arb. units)')
ax.set_ylabel('Relative Concentration')

# Add textbox with decay parameter and initial concentration
textstr = '\n'.join((
    # r'$\tau=%.2f$' % (par[1], ),
    r'$x_0={}$'.format([str(x) for x in init_conc], ),
    r'Rate: {}'.format([str(x) for x in rates], )))

props = dict(boxstyle='round', facecolor=colors[0], alpha=0.5)
ax.text(0.30, 0.10, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax.legend()
ax.grid(color="#424656", alpha=0.1)

plt.title("Example of the full solution to the rate equations")
# plt.subplots_adjust(top=0.93, bottom=0.09, left=0.06, right=0.98)
plt.show()

