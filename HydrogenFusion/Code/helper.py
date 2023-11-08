"""
A source file containing helper functions that can be used
across the project.
"""
import numpy as np
from scipy.optimize import curve_fit

def exp(x, a, b):
    """
    Exponential function.
    """
    return a * np.exp(b * x)


def quintic(x, a, b, c, d, e, f):
    """
    Quintic function.
    """
    return a + b * x + c * x**2 + d * x**3 + e * x**4 + f * x**5


def fit_func(x, y, func):
    """
    Fit a function to data.
    """
    popt, pcov = curve_fit(func, x, y)
    curve = func(x, *popt)
    return curve, popt, pcov


def find_decay_parameter(tau, sol, return_curve_fit=False):
    """
    Find the decay parameter of a solution by
    fitting an exponential function to it. This will
    help us quantify the rate of change of the solution.
    """
    # Fit an exponential function to the solution
    popt, pcov = curve_fit(exp, tau, sol)

    # Return the decay parameter
    if return_curve_fit:
        curve = exp(tau, *popt)
        return curve, popt, pcov
    else:
        return popt[1]
    

def distance(a, b):
    # Make weighted sum where points closer to the origin are weighted more
    # heavily
    return np.sum(np.abs(a - b) / np.abs(b))

def find_best_fit(func, data):
    # Define the range of parameters to search over
    k_range = np.linspace(0, 1000, 10)
    m_range = np.linspace(0, 1000, 10)

    # Initialize the best fit parameters
    best_k = 0
    best_m = 0

    # Initialize the best fit distance
    best_dist = np.inf

    # Loop over all parameters
    for i, k in enumerate(k_range):
        for j, m in enumerate(m_range):
            print(f"Now calculating {i}/{len(k_range)} and {j}/{len(m_range)}")
            # Calculate the fit
            fit = func(k, m)

            # Calculate the distance between the fit and the data
            dist = distance(fit, data)

            # Update the best fit parameters
            if dist < best_dist:
                best_dist = dist
                best_k = k
                best_m = m

    return best_k, best_m
