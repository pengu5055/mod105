"""
Fit the given arbitrary function to the data.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmasher as cmr
import palettable as pl
from scipy.optimize import curve_fit
from helper import *
from hydrogenfusion import HydrogenFusion

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


# Normalize the solutions to the sum of the initial concentrations
sol.y = sol.y / np.sum(init_conc)

H2 = sol.y[0]
Br2 = sol.y[1]
HBr = sol.y[2]
H = sol.y[3]
Br = sol.y[4]

derivative_HBr = np.gradient(HBr, t_eval)

RHS = lambda m: (H2 * Br2**(3/2)) / (m*Br2 + HBr)
arb_fit = lambda k, m: k * (H2 * Br2**(3/2)) / (m*Br2 + HBr)
scipy_fit = lambda t, k, m: k * (H2 * Br2**(3/2)) / (m*Br2 + HBr)

opt_k, opt_m = find_best_fit(arb_fit, derivative_HBr)

print(f"Optimal k: {opt_k}")
print(f"Optimal m: {opt_m}")

popt, pcov = curve_fit(scipy_fit, t_eval, derivative_HBr)
fit = scipy_fit(t_eval, *popt)


# Add textbox with ghetto fit parameters
textstr = '\n'.join((
    r'Ghetto Fit Parameters',
    r'$k=%.4e$' % (opt_k, ),
    r'$m=%.4e$' % (opt_m, )))
props = dict(boxstyle='round', facecolor='#eba646', alpha=0.5)
plt.text(0.98, 0.70, textstr, transform=plt.gca().transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment="right", bbox=props)
textstr = '\n'.join((
    r'SciPy Fit Parameters',
    r'$k=%.4e$' % (popt[0], ),
    r'$m=%.4e$' % (popt[1], )))
props = dict(boxstyle='round', facecolor='#9750a1', alpha=0.5)
plt.text(0.98, 0.53, textstr, transform=plt.gca().transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment="right", bbox=props)

k_hand = 1
m_hand = 0.44

textstr = '\n'.join((
    r'Hand Fit Parameters',
    r'$k=%.4e$' % (k_hand, ),
    r'$m=%.4e$' % (m_hand, )))
props = dict(boxstyle='round', facecolor='#649fe3', alpha=0.5)
plt.text(0.98, 0.31 + 0.05, textstr, transform=plt.gca().transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment="right", bbox=props)

plt.plot(t_eval, arb_fit(opt_k, opt_m), color='#eba646', label="Ghetto Fit")
plt.plot(t_eval, fit, color='#9750a1', label="SciPy Fit")
plt.plot(t_eval, arb_fit(k_hand, m_hand), color='#649fe3', label="Hand Fit")
plt.plot(t_eval, derivative_HBr, ls='--', color='#e83177', label="Data")
plt.xlabel("Time [arb. units]")
plt.ylabel("Derivative of $HBr$ concentration")
plt.legend()
plt.title("$HBr$ Concentration rate: Hand Fit vs. Ghetto Fit vs. SciPy")
plt.ylim(-0.02, 0.2)
plt.show()