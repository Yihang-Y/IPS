import numpy as np
from typing import *
import matplotlib.pyplot as plt

def leapfrog_integrate(x0: np.ndarray, p0: np.ndarray, n_steps: int, h: float, f: callable) -> list:
    """JUST ONE PARTICLE. Integrate systems of the form x'' = f(x), i.e. x' = p, p' = f(x).

    Args:
        x0 (np.ndarray): _description_
        p0 (np.ndarray): _description_
        n_steps (int): _description_
        h (float): _description_
        f (Callable): _description_

    Returns:
        List[np.ndarray]: x: the position trajectory of the point
    """
    d = len(x0)
    x = np.zeros((n_steps + 1, d), dtype = float)
    p = np.zeros((n_steps + 1, d), dtype = float)

    x[0], p[0] = x0, p0
    for i in range(n_steps):
        p_half = p[i] + h/2 * f(x[i])
        x[i+1] = x[i] + h * p_half
        p[i+1] = p_half + h/2 * f(x[i+1])

    return [x]

def gone_unstable(trajectories: list) -> bool:
    """Checks if any of the particles in the trajectories has shot off to infinity"""
    for particle_traject in trajectories:
        # go to the last point in the trajectory and see if it has flown off to infty
        if np.any(np.isnan(particle_traject[-1])):
            return True
        
    return False

def get_stability_threshold(integrator: callable, h_init: float = 1.0, n_integration_steps: int = 10000, max_iter: int = 20, verbose: bool = False) -> float:
    """_summary_

    Args:
        integrator (callable): Must be a function integrator(h, n_steps) that returns a list of position trajectories, each in the shape (n_steps, D) where D is the dimension.
        h_init (float, optional): _description_. Defaults to 1.0.
        n_integration_steps (int, optional): _description_. Defaults to 10000.
        max_iter (int, optional): _description_. Defaults to 20.
        verbose (bool, optional): _description_. Defaults to False.

    Returns:
        float: The threshold stepsize for stability
    """
    h = h_init
    n = 0

    # step 1: run simulation. Either it's stable or it's not
    if verbose: print("Step 1: big jumps")
    trajects = integrator(h, n_integration_steps)
    if gone_unstable(trajects):
        # make step size go downwards
        while gone_unstable(trajects) and n < max_iter:
            h /= 2
            if verbose: print(f"Unstable, trying {h=}")
            trajects = integrator(h, n_integration_steps)
            n += 1

        h_l, h_r = h, 2*h
    else:
        # make step size go upwards
        while not gone_unstable(trajects) and n < max_iter:
            h *= 2
            if verbose: print(f"Stable, trying {h=}")
            trajects = integrator(h, n_integration_steps)
            n += 1

        h_l, h_r = h/2, h

    if n == max_iter:
        if h > h_init:
            raise ValueError(f"Method was stable for all step sizes from {h_init} up to {h:.2e}")
        else:
            raise ValueError(f"Method was unstable for all step sizes from {h_init} down to {h:.2e}")

    if verbose:
        print(f"Step 1 complete, {h_l=}, {h_r=}")        

    # step 2: zero in on a better estimate of critical h. by this stage we know that h_l is stable and h_r is unstable
    n = 0 # reset our max iter counter
    while n < max_iter: # using a while loop so i can maybe add a different stopping criterion later
        h_mid = (h_l + h_r) / 2
        trajects = integrator(h_mid, n_integration_steps)
        if gone_unstable(trajects):
            if verbose: print(f"{h_mid=}, unstable")
            h_r = h_mid
        else:
            if verbose: print(f"{h_mid=}, stable")
            h_l = h_mid
            
        n += 1
            
    return (h_l + h_r) / 2

