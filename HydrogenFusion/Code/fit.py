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

