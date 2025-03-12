import numpy as np
from typing import *
import matplotlib.pyplot as plt

def leapfrog_integrate(x0: np.ndarray, p0: np.ndarray, n_steps: int, h: float, f: callable) -> tuple[np.ndarray, np.ndarray]:
    """Integrate systems of the form x'' = f(x), i.e. x' = p, p' = f(x).

    Args:
        x0 (np.ndarray): _description_
        p0 (np.ndarray): _description_
        n_steps (int): _description_
        h (float): _description_
        f (Callable): _description_

    Returns:
        Tuple[np.ndarray, np.ndarray]: x, p: the position and velocity trajectories
    """
    d = len(x0)
    x = np.zeros((n_steps + 1, d), dtype = float)
    p = np.zeros((n_steps + 1, d), dtype = float)

    x[0], p[0] = x0, p0
    for i in range(n_steps):
        p_half = p[i] + h/2 * f(x[i])
        x[i+1] = x[i] + h * p_half
        p[i+1] = p_half + h/2 * f(x[i+1])

    return x, p

def get_stability_threshold(f: callable):
    # this should be largely independent of the initial conditions, so we will always let these be the same
    h = 1
    # idea: keep doubling/halving h until I see the stability change. Then, do a kind of binary search to narrow down the range for h
