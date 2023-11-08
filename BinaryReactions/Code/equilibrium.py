"""
Exact solution for the three parameters specified
in the given task.
"""
"""
Basic test of the BinaryReactions class.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple
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
custom_colors = ["#845ec2", "#c031b5", "#cf3d2a"]
# A lighter version of custom_colors
box_colors = ["#a8a8c0", "#d17fb0", "#d88a7a"]

# Rather plot B + C

l_A1, = ax[0].plot(tau, sol.y[0], c=custom_colors[0], ls='--')
l_Astar1, = ax[0].plot(tau, sol.y[1], c=custom_colors[0], ls=':')
l_BC1, = ax[0].plot(tau, sol.y[2] + sol.y[3], c=custom_colors[0])

l_A2, = ax[0].plot(tau, sol2.y[0], c=custom_colors[1], ls='--')
l_Astar2, = ax[0].plot(tau, sol2.y[1], c=custom_colors[1], ls=':')
l_BC2, = ax[0].plot(tau, sol2.y[2]+sol2.y[3], c=custom_colors[1])

l_A3, = ax[0].plot(tau, sol3.y[0], c=custom_colors[2], ls='--')
l_Astar3, = ax[0].plot(tau, sol3.y[1], c=custom_colors[2], ls=':')
l_BC3, = ax[0].plot(tau, sol3.y[2]+sol3.y[3], c=custom_colors[2])

# Add legend
lgnd = ax[0].legend([(l_A1, l_A2, l_A3),
                    (l_Astar1, l_Astar2, l_Astar3),
                    (l_BC1, l_BC2, l_BC3),],
                    [r'$A$', r'$A^*$', r'$B+C$'],
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc='upper left', ncols=1, fontsize=10)

ax[0].grid(color="#424656", alpha=0.1)

# Add textboxes with parameters
textstr = r'$s_1=%.2f$' % (s, )
props = dict(boxstyle='round', facecolor=custom_colors[0], alpha=0.5)
ax[0].text(0.80, 0.60, textstr, transform=ax[0].transAxes, fontsize=11,
        verticalalignment='top', bbox=props)
textstr = r'$s_2=%.2f$' % (s2, )
props = dict(boxstyle='round', facecolor=custom_colors[1], alpha=0.5)
ax[0].text(0.80, 0.50, textstr, transform=ax[0].transAxes, fontsize=11,
        verticalalignment='top', bbox=props)
textstr = r'$s_3=%.2f$' % (s3, )
props = dict(boxstyle='round', facecolor=custom_colors[2], alpha=0.5)
ax[0].text(0.965, 0.40, textstr, transform=ax[0].transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment='right', bbox=props)

ax[0].set_ylabel('Relative Concentration')
ax[0].set_xlabel('Time (arb. units)')


# --- Continuous plot ---
solutions = []

# Define the parameters s
s = np.linspace(0.01, 1, 100)  # NOTE: After debug purposes, change to 100

# Solve the system for each s
for i in range(len(s)):
    print(f"Calculating solution {i+1}/{len(s)}")
    br = BinaryReaction(init_conc, tau, k, s[i])
    solutions.append(br.solve())

# Plot the results.
colors = cmr.take_cmap_colors(pl.scientific.sequential.Acton_10.mpl_colormap.reversed(),
                                len(solutions), cmap_range=(0., 0.7), return_fmt='hex')

for i in range(len(solutions)):
    ax[1].plot(tau, solutions[i].y[2] + solutions[i].y[3], c=colors[i])

ax[1].grid(color="#424656", alpha=0.1)
ax[1].set_ylabel('Relative Concentration')
ax[1].set_xlabel('Time (arb. units)')

norm = mpl.colors.LogNorm(vmin=s[0], vmax=s[-1])
sm = plt.cm.ScalarMappable(cmap=pl.scientific.sequential.Acton_10.mpl_colormap, norm=norm)

cbar = fig.colorbar(sm, ax=ax[1], orientation='vertical', pad=0.01)
cbar.set_label(r'$s$')



plt.suptitle("Solutions of the Binary reaction for different values of $s$")
plt.subplots_adjust(left=0.06, right=0.98, bottom=0.10, top=0.92)
plt.show()
