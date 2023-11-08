"""
Plot all three HydrogenFusion models on the same plot.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple
import palettable as pl
import cmasher as cmr
from hydrogenfusion import HydrogenFusion, HydrogenFusionStationaryBr, HydrogenFusionStationaryH

# Define the initial concentrations.
init_conc = np.array([1, 1, 0, 0, 0])

# Define the time points.
tau = np.linspace(0, 1, 1000)

# Define the constants.
rates = np.array([10, 1, 10, 0.1, 20])

# Initialize the HydrogenFusion objects.
hf = HydrogenFusion(init_conc, tau, rates)
hf_H = HydrogenFusionStationaryH(init_conc, tau, rates)

# Solve the system.
sol = hf.solve()
sol_H = hf_H.solve()

# Plot the results.
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Define the colors.
custom_colors = ["#9750a1", "#00c9bf", "#48817c"]
cm = pl.colorbrewer.sequential.Blues_7.mpl_colormap
colors = cmr.take_cmap_colors(cm, 6, cmap_range=(0.2, 1), return_fmt='hex')

cm_H = pl.colorbrewer.sequential.BuGn_7.mpl_colormap
colors_H = cmr.take_cmap_colors(cm_H, 6, cmap_range=(0.2, 1), return_fmt='hex')

# Normalize all solutions.
sol.y = sol.y / np.sum(init_conc)
sol_H.y = sol_H.y / np.sum(init_conc)

# Plot the solutions.
mod1_H2 = ax[0].plot(sol.t, sol.y[0], color=colors[0], linestyle='--', label=r'$\mathrm{H}_2$')
mod1_Br2 = ax[0].plot(sol.t, sol.y[1], color=colors[1], linestyle=':', label=r'$\mathrm{Br}_2$')
mod1_HBr = ax[0].plot(sol.t, sol.y[2], color=colors[2], label=r'$\mathrm{HBr}$')
mod1_H = ax[0].plot(sol.t, sol.y[3], color=colors[3], linestyle=(0, (5, 10)), label=r'$\mathrm{H}$')
mod1_Br = ax[0].plot(sol.t, sol.y[4], color=colors[4], linestyle=(0, (3, 10, 1, 10, 1, 10)), label=r'$\mathrm{Br}$')
mod1_sum = ax[0].plot(sol.t, np.sum(sol.y, axis=0), color=colors[5], linestyle='-', lw=1, label=r'$\sum_i c_i$')

mod2_H2 = ax[1].plot(sol_H.t, sol_H.y[0], color=colors_H[0], linestyle='--', label=r'$\mathrm{H}_2$')
mod2_Br2 = ax[1].plot(sol_H.t, sol_H.y[1], color=colors_H[1], linestyle=':', label=r'$\mathrm{Br}_2$')
mod2_HBr = ax[1].plot(sol_H.t, sol_H.y[2], color=colors_H[2], label=r'$\mathrm{HBr}$')
mod2_H = ax[1].plot(sol_H.t, sol_H.y[3], color=colors_H[3], linestyle=(0, (5, 10)), label=r'$\mathrm{H}$')
mod2_Br = ax[1].plot(sol_H.t, sol_H.y[4], color=colors_H[4], linestyle=(0, (3, 10, 1, 10, 1, 10)), label=r'$\mathrm{Br}$')
mod2_sum = ax[1].plot(sol_H.t, np.sum(sol_H.y, axis=0), color=colors_H[5], linestyle='-', lw=1, label=r'$\sum_i c_i$')

for x in ax:
    x.set_xlabel("Time [arb. units]")
    x.set_ylabel("Relative concentration")
    x.grid(color='#424656', alpha=0.1)
    x.hlines(1, 0, 1, color='#424656', linestyle='--')

ax[0].legend(loc='upper left', ncol=2, frameon=True)
ax[1].legend(loc='upper left', ncol=2, frameon=True)
ax[0].set_title("Exact solution")
ax[1].set_title(r"Stationary solution for $\dot{c}_{\mathrm{H}} = 0$")

# Add textboxes with the initial concentrations and rates.
textstr = '\n'.join((
    # r'$\tau=%.2f$' % (par[1], ),
    r'$x_0={}$'.format([str(x) for x in init_conc], ),
    r'Rate: {}'.format([str(x) for x in rates], )))

props = dict(boxstyle='round', facecolor=colors[0], alpha=0.5)
ax[0].text(0.98, 0.70, textstr, transform=ax[0].transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment="right", bbox=props)

textstr = '\n'.join((
    # r'$\tau=%.2f$' % (par[1], ),
    r'$x_0={}$'.format([str(x) for x in init_conc], ),
    r'Rate: {}'.format([str(x) for x in rates], )))

props = dict(boxstyle='round', facecolor=colors_H[0], alpha=0.5)
ax[1].text(0.98, 0.70, textstr, transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment="right", bbox=props)

plt.subplots_adjust(left=0.06, right=0.98, bottom=0.1, top=0.9)
plt.show()
