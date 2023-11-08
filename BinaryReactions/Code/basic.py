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
tau = np.linspace(0, 20, 1000)

# Define the constants.
k = 1000
s = 0.1

# Initialize the BinaryReaction object.
br = BinaryReaction(init_conc, tau, k, s)

# Solve the system.
sol = br.solve()

k2 = 1
s2 = 0.1

br2 = BinaryReaction(init_conc, tau, k2, s2)
sol2 = br2.solve()


# Plot the results.
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
colors = pl.cartocolors.sequential.TealGrn_7.hex_colors

ax[0].plot(tau, sol.y[0], c=colors[2], label='$A$', ls='--')
ax[0].plot(tau, sol.y[1], c=colors[3], label='$A^*$', ls=':')
ax[0].plot(tau, sol.y[2], c=colors[4], label='$B$', ls='-.')
ax[0].plot(tau, sol.y[3], c=colors[5], label='$C$')
ax[0].legend()

# Add textbox with parameters
textstr = '\n'.join((
    r'$k=%.2f$' % (k, ),
    r'$s=%.2f$' % (s, )))
props = dict(boxstyle='round', facecolor=colors[1], alpha=0.5)
ax[0].text(0.40, 0.20, textstr, transform=ax[0].transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

ax[0].set_ylabel('Relative Concentration')
ax[0].set_xlabel('Time (arb. units)')

ax[1].plot(tau, sol2.y[0], c=colors[2], label='$A$', ls='--')
ax[1].plot(tau, sol2.y[1], c=colors[3], label='$A^*$', ls=':')
ax[1].plot(tau, sol2.y[2], c=colors[4], label='$B$', ls='-.')
ax[1].plot(tau, sol2.y[3], c=colors[5], label='$C$')
ax[1].legend()

# Add textbox with parameters
textstr = '\n'.join((
    r'$k=%.2f$' % (k2, ),
    r'$s=%.2f$' % (s2, )))
props = dict(boxstyle='round', facecolor=colors[1], alpha=0.5)
ax[1].text(0.45, 0.20, textstr, transform=ax[1].transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

ax[1].set_ylabel('Relative Concentration')
ax[1].set_xlabel('Time (arb. units)')

plt.suptitle("Example Binary Reaction solutions")
plt.subplots_adjust(left=0.06, right=0.98, bottom=0.10, top=0.92)
plt.show()
