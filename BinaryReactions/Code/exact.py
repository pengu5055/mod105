"""
Exact solution for the three parameters specified
in the given task.
"""
"""
Basic test of the BinaryReactions class.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import palettable as pl
import cmasher as cmr
from binaryreaction import BinaryReaction


# Define the initial concentrations.
init_conc = np.array([1, 0, 0, 0])

# Define the time points.
tau = np.linspace(0, 20, 1000)

# Define the constants.
k = 1000
s = 0.1

# Initialize the BinaryReaction object.
br = BinaryReaction(init_conc, tau, k, s)

# Solve the system.
sol = br.solve()

# Repeat for other two parameters
s2 = 1
s3 = 10
br = BinaryReaction(init_conc, tau, k, s2)
sol2 = br.solve()
br = BinaryReaction(init_conc, tau, k, s3)
sol3 = br.solve()

# NOTE: Idea, solve for a range and do a continuous
# subfigure plot next to the main plot.

# Plot the results.
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
colors_1 = pl.cartocolors.sequential.TealGrn_7.hex_colors
colors_2 = cmr.take_cmap_colors(pl.cartocolors.sequential.PurpOr_7.mpl_colormap,
                                7, cmap_range=(0.2, 1.0), return_fmt='hex')
colors_3 = cmr.take_cmap_colors(pl.colorbrewer.sequential.YlOrRd_7.mpl_colormap,
                                7, cmap_range=(0.2, 1.0), return_fmt='hex')

l_A1 = ax[0].plot(tau, sol.y[0], c=colors_1[2], ls='--')
ax[0].plot(tau, sol.y[1], c=colors_1[3], ls=':')
ax[0].plot(tau, sol.y[2], c=colors_1[4], ls='-.')
ax[0].plot(tau, sol.y[3], c=colors_1[5])

l_A2 = ax[0].plot(tau, sol2.y[0], c=colors_2[2], ls='--')
ax[0].plot(tau, sol2.y[1], c=colors_2[3], ls=':')
ax[0].plot(tau, sol2.y[2], c=colors_2[4], ls='-.')
ax[0].plot(tau, sol2.y[3], c=colors_2[5])

l_A3 = ax[0].plot(tau, sol3.y[0], c=colors_3[2], ls='--')
ax[0].plot(tau, sol3.y[1], c=colors_3[3], ls=':')
ax[0].plot(tau, sol3.y[2], c=colors_3[4], ls='-.')
ax[0].plot(tau, sol3.y[3], c=colors_3[5])


ax[0].legend()
ax[0].grid(color="#424656", alpha=0.1)

# Add textbox with parameters
textstr = '\n'.join((
    r'$k=%.2f$' % (k, ),
    r'$s=%.2f$' % (s, )))
props = dict(boxstyle='round', facecolor=colors_1[1], alpha=0.5)
ax[0].text(0.40, 0.20, textstr, transform=ax[0].transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

ax[0].set_ylabel('Relative Concentration')
ax[0].set_xlabel('Time (arb. units)')


plt.suptitle("Example Binary Reaction solutions")
plt.subplots_adjust(left=0.06, right=0.98, bottom=0.10, top=0.92)
plt.show()
