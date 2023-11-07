"""
Basic evaluation of the Hydrogen Fusion process.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmasher as cmr
import palettable as pl
from hydrogenfusion import HydrogenFusion
from helper import *

# Define the initial concentrations
init_conc = np.array([1, 1, 0, 0, 0])
init_conc2 = np.array([4, 1, 0, 0, 3])

# Define the rate coefficients
rates = np.array([1, 1, 1, 1, 1])
rates2 = np.array([0.01, 2, 1, 0.1, 1])

# Define the time points to evaluate the rate equations
t_eval = np.linspace(0, 10, 1000)

# Initialize the HydrogenFusion class
hf = HydrogenFusion(init_conc, t_eval, rates)

# Solve the rate equations
sol = hf.solve()

# Repeat for the second set of initial concentrations
hf2 = HydrogenFusion(init_conc2, t_eval, rates2)
sol2 = hf2.solve()


# Plot the results
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

cm = pl.colorbrewer.sequential.BuPu_7.mpl_colormap
colors = cmr.take_cmap_colors(cm, 7, cmap_range=(0.2, 1), return_fmt='hex')

# --- Plot 1 ---
ax[0].plot(sol.t, sol.y[0], label='$H_2$', color=colors[1], linestyle='--')
ax[0].plot(sol.t, sol.y[1], label='$Br_2$', color=colors[5], linestyle='--')
ax[0].plot(sol.t, sol.y[2], label='$HBr$', color=colors[3])
ax[0].plot(sol.t, sol.y[3], label='$H$', color=colors[2], linestyle=':')
ax[0].plot(sol.t, sol.y[4], label='$Br$', color=colors[4], linestyle=':')

# NOTE: An idea is to make the curves dimensionless by normalizing to the sum
# of the initial concentrations (thus giving relative concentrations).

ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Concentration (M)')

# Add textbox with decay parameter and initial concentration
textstr = '\n'.join((
    # r'$\tau=%.2f$' % (par[1], ),
    r'$x_0={}$'.format([str(x) for x in init_conc], ),
    r'Rate: {}'.format([str(x) for x in rates], )))

props = dict(boxstyle='round', facecolor=colors[0], alpha=0.5)
ax[0].text(0.30, 0.10, textstr, transform=ax[0].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax[0].legend()
ax[0].grid(color="#424656", alpha=0.1)

# --- Plot 2 ---
ax[1].plot(sol2.t, sol2.y[0], label='$H_2$', color=colors[1], linestyle='--')
ax[1].plot(sol2.t, sol2.y[1], label='$Br_2$', color=colors[5], linestyle='--')
ax[1].plot(sol2.t, sol2.y[2], label='$HBr$', color=colors[3])
ax[1].plot(sol2.t, sol2.y[3], label='$H$', color=colors[2], linestyle=':')
ax[1].plot(sol2.t, sol2.y[4], label='$Br$', color=colors[4], linestyle=':')

ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('Concentration (M)')

# Add textbox with decay parameter and initial concentration
textstr = '\n'.join((
    # r'$\tau=%.2f$' % (par[1], ),
    r'$x_0={}$'.format([str(x) for x in init_conc2], ),
    r'Rate: {}'.format([str(x) for x in rates2], )))
props = dict(boxstyle='round', facecolor=colors[0], alpha=0.5)
ax[1].text(0.38, 0.22, textstr, transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax[1].legend()
ax[1].grid(color="#424656", alpha=0.1)

plt.suptitle("Example of the full solution to the rate equations")
plt.subplots_adjust(top=0.93, bottom=0.09, left=0.06, right=0.98)
plt.show()

