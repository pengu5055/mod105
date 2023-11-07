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