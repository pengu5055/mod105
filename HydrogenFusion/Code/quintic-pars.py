"""
Plot the dependence of the quintic fit parameters on the rate constants.
Probably best to make 2D heatmap for each parameter.
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

solutions = []
H2_range = np.linspace(0.1, 10, 100)
Br2_range = np.linspace(0.1, 10, 100)

for i, H2 in enumerate(H2_range):
    print(f"Progress: {i}/{len(H2_range)}")
    column = []
    for Br2 in Br2_range:
        rates = np.array([H2, Br2, 1, 1, 1])
        hf = HydrogenFusion(init_conc, t_eval, rates)
        sol = hf.solve()
        sol.y = sol.y / np.sum(init_conc)
        column.append(fit_func(sol.t, sol.y[2], quintic)[1][2])

    solutions.append(column)

solutions = np.array(solutions)
df = pd.DataFrame(solutions, index=H2_range, columns=Br2_range)
df.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='Fit', mode='w', complevel=9)

# Plot
fig, ax = plt.subplots(2, 3, figsize=(10, 10))

# Define the colors.
cm_a = pl.colorbrewer.sequential.BuPu_7.mpl_colormap
cm_b = pl.colorbrewer.sequential.BuGn_7.mpl_colormap
cm_c = pl.colorbrewer.sequential.OrRd_7.mpl_colormap
cm_d = pl.colorbrewer.sequential.PuBu_7.mpl_colormap
cm_e = pl.colorbrewer.sequential.PuRd_7.mpl_colormap

# DEBUG plot first subplot
ax[0].contourf(H2_range, Br2_range, solutions[1][0], levels=100, cmap=cm_a)


ax[0].set_xlabel(r'Initial concentration of $H_2$')
ax[0].set_ylabel(r'Initial concentration of $Br_2$')

plt.show()