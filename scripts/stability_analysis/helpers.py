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

def get_stability_threshold(f: callable, max_iter: int = 20, verbose: bool = False):
    # this should be largely independent of the initial conditions, so we will always let these be the same
    h = 1
    x0 = np.array([1])
    p0 = np.array([-1])
    n_steps = 10000
    n = 0

    # step 1: run simulation. Either it's stable or it's not
    x, p = leapfrog_integrate(x0, p0, n_steps, h, f)
    if np.isnan(x[-1]):
        # make step size go downwards
        while np.isnan(x[-1]):
            h /= 2
            x,p = leapfrog_integrate(x0, p0, n_steps, h, f)

        h_l, h_r = h, 2*h
    else:
        # make step size go upwards
        while not np.isnan(x[-1]):
            h *= 2
            x,p = leapfrog_integrate(x0, p0, n_steps, h, f)

        h_l, h_r = h/2, h

    if verbose:
        print(f"Step 1 complete, {h_l=}, {h_r=}")        

    # step 2: zero in on a better estimate of critical h. by this stage we know that h_l is stable and h_r is unstable
    while n < max_iter: # using a while loop so i can maybe add a different stopping criterion later
        h_mid = (h_l + h_r) / 2
        x,p = leapfrog_integrate(x0, p0, n_steps, h_mid, f)
        if np.isnan(x[-1]):
            if verbose: print(f"{h_mid=}, unstable")
            h_r = h_mid
        else:
            if verbose: print(f"{h_mid=}, stable")
            h_l = h_mid
            
        n += 1
            
    return (h_l + h_r) / 2

