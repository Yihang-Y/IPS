# NOT WORKING

import numpy as np
import scipy as sp
import sys
sys.path.append('../../build')
sys.path.append('../')
import IPSModule as ips
import multiprocessing as mp

eps = 1.0
sigma = 1.0
rad = 5.0
T = 0.1
gamma = 1.0
n_simulations = 10000
n_steps = 10000
dt = 0.001
bins = 10

# set up the pair force and confinement
pair_force_config = {
    "type": "LennardJones",
    "eps": eps,
    "sigma": sigma
}
confinement_config = {
    "type": "Radial",
    "rad": rad
}


def U_LJ(q1, q2):
    r = np.linalg.norm(q1-q2, 2)
    return 4*eps*((sigma/r)**12 - (sigma/r)**6)


def U_confine(q1, q2):
    u_conf = 0
    
    # loop over each particle
    for q in [q1, q2]:
        # check if particle is in the naughty zone
        if np.linalg.norm(q) >= rad:
            # loop over each dimension
            for i in range(len(q1)):
                u_conf += (-rad-1)*np.log(np.abs(q[i] - rad - 1) - (q[i]**2 + 2*q[i])/2)

    return u_conf

def calc_rho(q1, q2, beta):
    # probability density
    if np.linalg.norm(q1, 2) > (rad + 1) or np.linalg.norm(q2, 2) > (rad + 1):
        # it is impossible for a particle starting inside the confinement radius to end up outside of it
        # is this true? Maybe because of noise in the SDE the expected density is nonzero outside the confinement zone?
        rho_val = 1e-16 # avoid divide by zero
    elif np.linalg.norm(q1-q2, 2) == 0:
        rho_val = 1e16 # essentially infinity, but we do not actually pass to U_LJ because then we'd get a divide by zero
    else:
        rho_val = np.exp(-beta*U_LJ(q1, q2) + U_confine(q1, q2))

    return rho_val

def run_simulation():
    """Returns 4D array"""
    particles = ips.LangevinSystem(2, gamma, T)
    particles.positions = np.random.uniform(-np.sqrt(rad), -np.sqrt(rad), size = (2,2)).tolist()
    simulator = ips.IPS_Simulator_Langevin(particles)
    simulator.init(pair_force_config, confinement_config)
    simulator.simulate_n_steps(dt, n_steps)
    q1, q2 = particles.get_positions()
    q = np.array([q1, q2]).flatten()

    return q

if __name__ == "__main__":
    # end_Q = np.zeros((n_simulations, 4))

    # # parallelise this
    # for n in n_simulations:
    #     end_Q[n] = run_simulation()

    # hist, _ = np.histogramdd(end_Q, bins = 10, range = [(-(rad+1),rad+1)]*4)
    
    rho = np.zeros((bins,bins,bins,bins))
    # now compare to the actual distribution
    qs = np.linspace(-(rad+1), rad+1, 10)
    Q1, Q2, Q3, Q4 = np.meshgrid(qs, qs, qs, qs)
    for i in range(len(qs)):
        for j in range(len(qs)):
            for k in range(len(qs)):
                for l in range(len(qs)):
                    q1 = np.array([qs[i], qs[j]])
                    q2 = np.array([qs[k], qs[l]])
                    rho[i,j,k,l] = calc_rho(q1, q2, beta = 1/T)

    entropy = sp.stats.entropy(rho.flatten(), rho.flatten())
    print(f"{entropy=}")