"""
Contains helper functions to be used everywhere.
"""
import numpy as np

def find_stop(sol: np.ndarray) -> int:
    """
    Find the point where the reaction stops.
    """
    for i, val in enumerate(sol):
        if np.abs(val) < 1e-4:
            return i
    return -1
