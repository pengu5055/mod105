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

# This solution is terrible, but I'm in a hurry.
solutions_a = []
solutions_b = []
solutions_c = []
solutions_d = []
solutions_e = []
solutions_f = []
H2_range = np.linspace(0.1, 10, 100)
Br2_range = np.linspace(0.1, 10, 100)

if True:
    for i, H2 in enumerate(H2_range):
        print(f"Progress: {i}/{len(H2_range)}")
        # Ditto.
        column_a = []
        column_b = []
        column_c = []
        column_d = []
        column_e = []
        column_f = []
        for Br2 in Br2_range:
            rates = np.array([H2, Br2, 1, 1, 1])
            hf = HydrogenFusion(init_conc, t_eval, rates)
            sol = hf.solve()
            sol.y = sol.y / np.sum(init_conc)
            column_a.append(fit_func(sol.t, sol.y[2], quintic)[1][0])
            column_b.append(fit_func(sol.t, sol.y[2], quintic)[1][1])
            column_c.append(fit_func(sol.t, sol.y[2], quintic)[1][2])
            column_d.append(fit_func(sol.t, sol.y[2], quintic)[1][3])
            column_e.append(fit_func(sol.t, sol.y[2], quintic)[1][4])
            column_f.append(fit_func(sol.t, sol.y[2], quintic)[1][5])
    
        solutions_a.append(column_a)
        solutions_b.append(column_b)
        solutions_c.append(column_c)
        solutions_d.append(column_d)
        solutions_e.append(column_e)
        solutions_f.append(column_f)
    
    solutions_a = np.array(solutions_a)
    df_a = pd.DataFrame(solutions_a, index=H2_range, columns=Br2_range)
    df_a.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='a', mode='w', complevel=9)
    df_b = pd.DataFrame(solutions_b, index=H2_range, columns=Br2_range)
    df_b.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='b', mode='a', complevel=9)
    df_c = pd.DataFrame(solutions_c, index=H2_range, columns=Br2_range)
    df_c.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='c', mode='a', complevel=9)
    df_d = pd.DataFrame(solutions_d, index=H2_range, columns=Br2_range)
    df_d.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='d', mode='a', complevel=9)
    df_e = pd.DataFrame(solutions_e, index=H2_range, columns=Br2_range)
    df_e.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='e', mode='a', complevel=9)
    df_f = pd.DataFrame(solutions_f, index=H2_range, columns=Br2_range)
    df_f.to_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='f', mode='a', complevel=9)
    quit()

# Load the data.
df = pd.read_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='a')
solutions_a = df.values
df = pd.read_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='b')
solutions_b = df.values
df = pd.read_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='c')
solutions_c = df.values
df = pd.read_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='d')
solutions_d = df.values
df = pd.read_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='e')
solutions_e = df.values
df = pd.read_hdf('./HydrogenFusion/Results/quintic-pars.h5', key='f')
solutions_f = df.values


# Plot
fig, ax = plt.subplots(2, 3, figsize=(10, 10))

# Define the colors.
cm_a = pl.colorbrewer.sequential.BuPu_7.mpl_colormap
cm_b = pl.colorbrewer.sequential.BuGn_7.mpl_colormap
cm_c = pl.colorbrewer.sequential.OrRd_7.mpl_colormap
cm_d = pl.colorbrewer.sequential.PuBu_7.mpl_colormap
cm_e = pl.colorbrewer.sequential.PuRd_7.mpl_colormap
cm_f = pl.colorbrewer.sequential.RdPu_7.mpl_colormap

# DEBUG plot first subplot
ax[0, 0].contourf(H2_range, Br2_range, solutions_a, levels=100, cmap=cm_a)


ax[0, 0].set_xlabel(r'Initial concentration of $H_2$')
ax[0, 0].set_ylabel(r'Initial concentration of $Br_2$')

plt.show()