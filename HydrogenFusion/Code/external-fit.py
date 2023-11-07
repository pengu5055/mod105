"""
Use an external program to fit a few arbitrary functions to the data
to see which one fits best. Basically must truncate 1000 data points
to only 20 significant points, or else I need to pay money.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import palettable as pl

# Load the data
df = pd.read_csv('./HydrogenFusion/Results/basic-sol.csv', index_col=0)

# Get the data
t = df.index.values
x = df['x'].values

# Truncate the data
t_t = t[::50]
x_t = x[::50]

# Plot the data
fig, ax = plt.subplots()

colors = pl.cartocolors.sequential.Purp_7.hex_colors

ax.plot(t, x, label='x', color=colors[6])
ax.scatter(t_t, x_t, color=colors[3], s=5, label="Truncated data")

ax.legend()
plt.show()

truncated_df = pd.DataFrame(x_t, columns=['x'], index=t_t).to_csv('./HydrogenFusion/Results/truncated-sol.csv')

# Got a = −0.03070366, b = 0.4031943, c = −0.02803624, d = −0.007001959, e = 0.00139959, f = −0.00006806583
# for a quintic regression. Coefficients are in the order of lowest to highest power.
