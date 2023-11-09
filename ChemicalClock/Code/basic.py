"""
Try to get a basic plot of oscillations in
the concentrations of Iodine and Iodide.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple
import palettable as pl
import cmasher as cmr
from chemicalclock import ChemicalClock

# Define function that searches for distance from 0.
# If distance is less than 1e-3, return index.

def find_stop(sol: np.ndarray) -> int:
    """
    Find the point where the reaction stops.
    """
    for i, val in enumerate(sol):
        if np.abs(val) < 1e-4:
            return i
    return -1

# Define the initial concentrations.
init_conc = np.array([1, 0, 10, 0, 0, 10, 0, 0])

# Define the time points.
tau = np.linspace(0, 40, 1000)

# Define the rate coefficients.
rate_factor = 1
fixed = 0.1
second_par = fixed/rate_factor
rates = np.array([rate_factor, second_par])
# p_fast_factor = 1
# q_fast_factor = 1

# Create the chemical clock.
clock = ChemicalClock(init_conc, tau, rates, p_fast_factor=100, q_fast_factor=100)

# Solve the model equations.
sol = clock.solve()

# Repeat for other two parameters
rate_factor2 = 10
second_par = fixed/rate_factor2
rates = np.array([rate_factor2, second_par])
clock = ChemicalClock(init_conc, tau, rates, p_fast_factor=100, q_fast_factor=100)
sol2 = clock.solve()

rate_factor3 = 100
second_par = fixed/rate_factor3
rates = np.array([rate_factor3, second_par])
clock = ChemicalClock(init_conc, tau, rates, p_fast_factor=100, q_fast_factor=100)
sol3 = clock.solve()

# Normalize the solutions.
sol.y = sol.y / np.max(sol.y, axis=1)[:, None]
sol2.y = sol2.y / np.max(sol2.y, axis=1)[:, None]
sol3.y = sol3.y / np.max(sol3.y, axis=1)[:, None]

# Plot the results.
fig, ax = plt.subplots(1, 1)

custom_colors = ["#845ec2", "#c031b5", "#cf3d2a"]

I1, = ax.plot(sol.t, sol.y[0], label=r'$[\mathrm{I^-}]$', color=custom_colors[0])
I2_1, = ax.plot(sol.t, sol.y[1], label=r'$[\mathrm{I_2}]$', color=custom_colors[0], linestyle='--')



I2, = ax.plot(sol.t, sol.y[2], label=r'$[\mathrm{I_3^-}]$', color=custom_colors[1])
I2_2, = ax.plot(sol.t, sol.y[3], label=r'$[\mathrm{I_5^-}]$', color=custom_colors[1], linestyle='--')

I3, = ax.plot(sol.t, sol.y[4], label=r'$[\mathrm{I_2O_5}]$', color=custom_colors[2])
I2_3, = ax.plot(sol.t, sol.y[5], label=r'$[\mathrm{I_2O_7}]$', color=custom_colors[2], linestyle='--')

# Find incidence point where [M] = [N]
# and plot the vertical line.
incidence = np.argmin(np.abs(sol.y[0] - sol.y[1]))
# ax.axvline(sol.t[incidence], color='k', linestyle='--', label='Equilibrium point', alpha=0.5)

stop = find_stop(np.gradient(sol.y[0]))
stop2 = find_stop(np.gradient(sol2.y[0]))
stop3 = find_stop(np.gradient(sol3.y[0]))
# ax.axvline(sol.t[stop], color='k', linestyle=':', label='Reaction stops', alpha=0.5)

ax.legend(((I1, I2, I3), (I2_1, I2_2, I2_3),),
          (r'$[\mathrm{I^-}]$', r'$[\mathrm{I_2}]$'),
          handler_map={tuple: HandlerTuple(ndivide=None)},
          loc='lower right', ncol=3, fontsize=11)

# Add textboxes with parameters
textstr = '\n'.join(
        (r'$\lambda_1=%.4e$' % (rate_factor, ),
         r'$\tau_1=%.4f$' % (stop, )))
props = dict(boxstyle='round', facecolor=custom_colors[0], alpha=0.5)
ax.text(0.98, 0.87, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment="right", bbox=props)
textstr = '\n'.join(
        (r'$\lambda_2=%.4e$' % (rate_factor2, ),
         r'$\tau_2=%.4f$' % (stop2, )))
props = dict(boxstyle='round', facecolor=custom_colors[1], alpha=0.5)
ax.text(0.98, 0.70, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment="right", bbox=props)
textstr = '\n'.join(
        (r'$\lambda_3=%.4e$' % (rate_factor3, ),
         r'$\tau_3=%.4f$' % (stop3, )))
props = dict(boxstyle='round', facecolor=custom_colors[2], alpha=0.5)
ax.text(0.98, 0.40, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment="right", bbox=props)


ax.set_xlabel('Time')
ax.set_ylabel('Relative Concentration')
ax.grid(color="#424656", alpha=0.1)
ax.set_title('Full solution for a Iodine Chemical Clock')
plt.show()
