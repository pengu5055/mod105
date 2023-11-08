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
tau = np.linspace(0, 10, 1000)

# Define the constants.
rates = np.array([1, 1, 1, 1, 1])

# Initialize the HydrogenFusion objects.
hf = HydrogenFusion(init_conc, tau, rates)
hf_Br = HydrogenFusionStationaryBr(init_conc, tau, rates)
hf_H = HydrogenFusionStationaryH(init_conc, tau, rates)

# Solve the system.
sol = hf.solve()
sol_Br = hf_Br.solve()
sol_H = hf_H.solve()

# Plot the results.
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Define the colors.
custom_colors = ["#9750a1", "#00c9bf", "#48817c"]
cm = pl.colorbrewer.sequential.Blues_7.mpl_colormap
colors = cmr.take_cmap_colors(cm, 6, cmap_range=(0.2, 1), return_fmt='hex')

cm_Br = pl.colorbrewer.sequential.RdPu_7.mpl_colormap
colors_Br = cmr.take_cmap_colors(cm_Br, 6, cmap_range=(0.2, 1), return_fmt='hex')

cm_H = pl.colorbrewer.sequential.BuGn_7.mpl_colormap
colors_H = cmr.take_cmap_colors(cm_H, 6, cmap_range=(0.2, 1), return_fmt='hex')

# --- Plot 1 ---

# Plot the solutions.
mod1_H2 = ax[0].plot(sol.t, sol.y[0], color=colors[0], linestyle='--')
mod1_Br2 = ax[0].plot(sol.t, sol.y[1], color=colors[1], linestyle=':')
mod1_HBr = ax[0].plot(sol.t, sol.y[2], color=colors[2])
mod1_H = ax[0].plot(sol.t, sol.y[3], color=colors[3], linestyle=(0, (5, 10)))
mod1_Br = ax[0].plot(sol.t, sol.y[4], color=colors[4], linestyle=(0, (3, 10, 1, 10, 1, 10)))
mod1_sum = ax[0].plot(sol.t, np.sum(sol.y, axis=0), color=colors[5], linestyle='-', lw=1)

mod2_H2 = ax[0].plot(sol_Br.t, sol_Br.y[0], color=colors_Br[0], linestyle='--')
mod2_Br2 = ax[0].plot(sol_Br.t, sol_Br.y[1], color=colors_Br[1], linestyle=':')
mod2_HBr = ax[0].plot(sol_Br.t, sol_Br.y[2], color=colors_Br[2])
mod2_H = ax[0].plot(sol_Br.t, sol_Br.y[3], color=colors_Br[3], linestyle=(0, (5, 10)))
mod2_Br = ax[0].plot(sol_Br.t, sol_Br.y[4], color=colors_Br[4], linestyle=(0, (3, 10, 1, 10, 1, 10)))
mod2_sum = ax[0].plot(sol_Br.t, np.sum(sol_Br.y, axis=0), color=colors_Br[5], linestyle='-', lw=1)

mod3_H2 = ax[0].plot(sol_H.t, sol_H.y[0], color=colors_H[0], linestyle='--')
mod3_Br2 = ax[0].plot(sol_H.t, sol_H.y[1], color=colors_H[1], linestyle=':')
mod3_HBr = ax[0].plot(sol_H.t, sol_H.y[2], color=colors_H[2])
mod3_H = ax[0].plot(sol_H.t, sol_H.y[3], color=colors_H[3], linestyle=(0, (5, 10)))
mod3_Br = ax[0].plot(sol_H.t, sol_H.y[4], color=colors_H[4], linestyle=(0, (3, 10, 1, 10, 1, 10)))
mod3_sum = ax[0].plot(sol_H.t, np.sum(sol_H.y, axis=0), color=colors_H[5], linestyle='-', lw=1)

ax[0].hlines(1, 0, 10, color='#424656', linestyle='--')

# --- Plot 2 ---

ax[1].plot(sol.t, sol.y[0], color=colors[0], linestyle='--')

plt.show()
