"""
Fit quintic polynomial to the data.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmasher as cmr
import palettable as pl
from helper import *
from hydrogenfusion import HydrogenFusion

# Do the usual initialization.
init_conc = np.array([1, 1, 0, 0, 0])
rates = np.array([1, 1, 1, 1, 1])
t_eval = np.linspace(0, 10, 1000)
hf = HydrogenFusion(init_conc, t_eval, rates)
sol = hf.solve()

# Normalize the solutions to the sum of the initial concentrations
sol.y = sol.y / np.sum(init_conc)

# Want to fit to the HBr concentration.
fit, par, err = fit_func(sol.t, sol.y[2], quintic)

# Plot
fig, ax = plt.subplots(1, 1, figsize=(5, 5))

# Define the colors.
cm = pl.colorbrewer.sequential.BuPu_7.mpl_colormap
colors = cmr.take_cmap_colors(cm, 7, cmap_range=(0.2, 1), return_fmt='hex')

# Truncate the data to say 20 significant points
t_trunc = sol.t[::50]
data_trunc = sol.y[2][::50]
ax.scatter(t_trunc, data_trunc, label='$HBr$', color=colors[3])
ax.plot(sol.t, fit, label='Quintic Fit', color=colors[5])

# Add textbox with quintic parameters
textstr = '\n'.join((
    r'$x^0:\quad a=%.4e$' % (par[0], ),
    r'$x^1:\quad b=%.4e$' % (par[1], ),
    r'$x^2:\quad c=%.4e$' % (par[2], ),
    r'$x^3:\quad d=%.4e$' % (par[3], ),
    r'$x^4:\quad e=%.4e$' % (par[4], ),
    r'$x^5:\quad f=%.4e$' % (par[5], )))
props = dict(boxstyle='round', facecolor=colors[0], alpha=0.5)
ax.text(0.43, 0.385, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

ax.set_xlabel('Time [arb. units]')
ax.set_ylabel('Relative Concentration')

ax.legend(loc='upper left', ncol=1, frameon=False)
ax.grid(color="#424656", alpha=0.1)
plt.title("Quintic Fit to $HBr$ Concentration")
plt.show()
