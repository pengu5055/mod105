"""
Plot the rate of synthesis of HBr as a function of time
for different initial concentrations of H2 and Br2.
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

init_conc1 = np.array([100, 1, 0, 0, 0])
init_conc2 = np.array([1, 1, 0, 0, 0])
init_conc3 = np.array([1, 100, 0, 0, 0])

hf1 = HydrogenFusion(init_conc1, t_eval, rates)
hf2 = HydrogenFusion(init_conc2, t_eval, rates)
hf3 = HydrogenFusion(init_conc3, t_eval, rates)

sol1 = hf1.solve()
sol2 = hf2.solve()
sol3 = hf3.solve()

# Normalize the solutions to the sum of the initial concentrations
sol1.y = sol1.y / np.sum(init_conc1)
sol2.y = sol2.y / np.sum(init_conc2)
sol3.y = sol3.y / np.sum(init_conc3)

# Plot
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# Define the colors.
custom_colors = ["#9750a1", "#00c9bf", "#48817c"]

ax[0].plot(sol1.t, sol1.y[2], color=custom_colors[0], label=r'$\frac{\mathrm{H}_2}{\mathrm{Br}_2} = 100$')
ax[0].plot(sol2.t, sol2.y[2], color=custom_colors[1], label=r'$\frac{\mathrm{H}_2}{\mathrm{Br}_2} = 1$')
ax[0].plot(sol3.t, sol3.y[2], color=custom_colors[2], label=r'$\frac{\mathrm{H}_2}{\mathrm{Br}_2} = 0.01$')

ax[0].set_xlabel('Time [arb. units]')
ax[0].set_ylabel('Relative Concentration')

ax[0].legend(loc='upper left', ncol=1, frameon=False)
ax[0].grid(color="#424656", alpha=0.1)

# --- Continuous plot ---
# Define the initial concentrations.
solutions = []
H2_range = np.linspace(0.001, 100, 100)
Br2_range = np.linspace(0.001, 100, 100)

for i, H2 in enumerate(H2_range):
    print(f"Progress: {i}/{len(H2_range)}")
    column = []
    for Br2 in Br2_range:
        init_conc = np.array([H2, Br2, 0, 0, 0])
        hf = HydrogenFusion(init_conc, t_eval, rates)
        sol = hf.solve()
        sol.y = sol.y / np.sum(init_conc)
        column.append(fit_func(sol.t, sol.y[2], quintic)[1][2])

    solutions.append(column)

solutions = np.array(solutions)
df = pd.DataFrame(solutions, index=H2_range, columns=Br2_range)
df.to_hdf('./HydrogenFusion/Results/synthesis.h5', key='Synthesis', mode='w', complevel=9)

# Plot
# Define the colors.
norm = mpl.colors.Normalize(vmin=0, vmax=1)
cm = pl.colorbrewer.sequential.BuPu_7.mpl_colormap

ax[1].contourf(H2_range, Br2_range, solutions, levels=100, cmap=cm, norm=norm)
ax[1].contour(H2_range, Br2_range, solutions, levels=20)

ax[1].set_xlabel(r'Initial concentration of $H_2$')
ax[1].set_ylabel(r'Initial concentration of $Br_2$')

sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
cbar = plt.colorbar(sm, ax=ax[1])
cbar.set_label('Size of the leading coef. in the quintic expansion for the rate of synthesis of $HBr$')

plt.show()


