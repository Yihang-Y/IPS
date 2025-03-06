import subprocess
from typing import *
import numpy as np

def parse_bmm_output(bmm_output: str) -> List:
    """Extract the trajectories from the BMM output."""
    trajectories = []
    current_traj = 0

    step_size_arr = []
    num_steps_arr = []

    lines = bmm_output.splitlines()

    for line in lines:
        line = line.strip()

        # we don't need these NUM_STEPS and STEP_SIZE right now because the caller in Python will already know that, but maybe in future this will be useful so we leave it in.
        if line.startswith("# NUM_STEPS"):
            num_steps = int(line.split()[-1])
            num_steps_arr.append(num_steps)
        elif line.startswith("# STEP_SIZE"):
            step_size = float(line.split()[-1])
            step_size_arr.append(step_size)
        else:
            values = list(map(float, line.split())) 
            trajectories.append(values)
            current_traj += 1

    return trajectories

def BrownianMotionMain(n_steps: int, step_size: float, n_trajects: int, path_to_exe: str, use_wsl: bool = False) -> List[str]:
    """Simple wrapper that calls C++ executable for 1D Brownian motion overdamped Langevin simulations, and returns results as a list of strings.

    Args:
        n_steps (int): number of steps
        step_size (float): step size
        n_trajects (int): number of trajectories
        path_to_exe (str): the path to the BMM executable
        use_wsl (bool, optional): whether or not to use WSL (needed for Windows, not for Linux (or mac?)). Defaults to False.

    Returns:
        List[str]: Result of the calculations (have a look at what BrownianMotionMain.exe outputs to see format).
    """

    try:
        if use_wsl:
            result = subprocess.run(["wsl", "-e", path_to_exe, str(n_steps), str(step_size), str(n_trajects)], shell = True, capture_output = True)
        else:
            result = subprocess.run([path_to_exe, str(n_steps), str(step_size), str(n_trajects)], capture_output = True)
    except Exception as e:
        raise Exception("Problem running the executable: ", e)
    
    if result.stderr:
        raise Exception("Error during calculations: ", result.stderr.decode("ascii"))
    else:
        decoded_output = result.stdout.decode("ascii")
        trajectories = parse_bmm_output(decoded_output)

        return trajectories
    

def read_ips_particle_data(output_dump, dim: int) -> List[np.array]:
    # read data from file
    # each line is the position of a particle at a given time
    # and between different time steps there is a blank line
    data = output_dump.split('\n\n')
        
    particle_data = []
    for block in data:
        block_data = []
        for line in block.split('\n'):
            try:
                arr = list(map(float, line.split(" ")))[:dim]
            except ValueError:
                print(f"Error reading line: {line}")
                continue
            block_data.append(arr)
        if block_data != []:
            particle_data.append(np.array(block_data))
    return particle_data

def IPS(n_steps: int, step_size: float, num_particles: int, dim: int, output_interval: int, path_to_exe: str, use_wsl: bool = False) -> List[str]:
    """Simple wrapper that calls C++ executable for the IPS, and returns results as a list of strings.

    Args:
        n_steps (int): number of steps
        step_size (float): step size
        num_particles (int): number of particles
        dim (int): The dimension of the particles (currently only d = 2 or 3 supported)
        output_interval (int): How often to output the particle info. So the simulation will still run for every step,
        but you will only get the trajectery info every x number of steps. Useful to save space and make plotting faster in Python.
        path_to_exe (str): the path to the IPS executable
        use_wsl (bool, optional): whether or not to use WSL (needed for Windows, not for Linux (or mac?)). Defaults to False.

    Returns:
        List[str]: Result of the calculations (have a look at what IPS.cpp outputs to see format).
    """

    try:
        if use_wsl:
            result = subprocess.run(["wsl", "-e", path_to_exe, str(n_steps), str(step_size), str(num_particles), str(dim), str(output_interval)], shell = True, capture_output = True)
        else:
            result = subprocess.run([path_to_exe, str(n_steps), str(step_size),  str(num_particles), str(dim), str(output_interval)], capture_output = True)
    except Exception as e:
        raise Exception("Problem running the executable: ", e)
    
    if result.stderr:
        raise Exception("Error during calculations: ", result.stderr.decode("ascii"))
    else:
        decoded_output = result.stdout.decode("ascii")
        trajectories = read_ips_particle_data(decoded_output, dim)

        return trajectories