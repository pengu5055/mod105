"""
Basic evaluation of the Hydrogen Fusion process.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import palettable as pl
from hydrogenfusion import HydrogenFusion
from helper import *

# Define the initial concentrations
init_conc = np.array([1, 1, 0, 0, 0])

# Define the rate coefficients
rates = np.array([1, 1, 1, 1, 1])

# Define the time points to evaluate the rate equations
t_eval = np.linspace(0, 10, 1000)

# Initialize the HydrogenFusion class
hf = HydrogenFusion(init_conc, t_eval, rates)

# Solve the rate equations
sol = hf.solve()

# pd.DataFrame(sol.y.T, columns=['u', 'v', 'x', 'y', 'z']).to_csv('./HydrogenFusion/Results/basic-sol.csv')
df = pd.DataFrame(sol.y[2], columns=['x'], index=sol.t).to_csv('./HydrogenFusion/Results/basic-sol.csv')



# Plot the results
fig, ax = plt.subplots()

colors = pl.cartocolors.sequential.Purp_7.hex_colors

fit, par, err = find_decay_parameter(sol.t, sol.y[2], return_curve_fit=True)

ax.plot(sol.t, fit, label='fit', color=colors[2])
ax.scatter(sol.t, sol.y[2], label='x', color=colors[6], s=4)

# NOTE: An idea is to make the curves dimensionless by normalizing to the sum
# of the initial concentrations (thus giving relative concentrations).

ax.set_xlabel('Time (s)')
ax.set_ylabel('Concentration (M)')
plt.title("Concentration of hydrogen bromide over time")

# Add textbox with decay parameter and initial concentration
textstr = '\n'.join((
    r'$\tau=%.2f$' % (par[1], ),
    r'$x_0={}$'.format(init_conc, )))

props = dict(boxstyle='round', facecolor=colors[0], alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

ax.legend()
plt.show()

