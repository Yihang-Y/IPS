# Cluster Formation Time Simulation (using IPSModule Langevin system)
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import IPSModule as ips

# -----------------------------
# Simulation parameters
# -----------------------------
gamma = 1.0
temperature = 0.5
dt = 0.0005
rad = 30.0
n_particles = 50
epsilon = 1.0
sigma = 1.0

# -----------------------------
# Generate random particle positions in circular region
# -----------------------------
def generate_random_particles(n_particles=50, box_radius=20.0):
    np.random.seed(None)
    while True:
        positions = np.random.uniform(-box_radius, box_radius, size=(n_particles, 2))
        distances = np.linalg.norm(positions, axis=1)
        if np.all(distances <= box_radius):
            return positions

# -----------------------------
# Count number of clusters using DBSCAN
# -----------------------------
def count_clusters(positions, eps=1.7, min_samples=3):
    if np.isnan(positions).any():
        print("Skipping frame due to NaN in positions")
        return 0  # Return 0 clusters if NaN detected
    db = DBSCAN(eps=eps, min_samples=min_samples)
    labels = db.fit_predict(positions)
    return len(set(labels)) - (1 if -1 in labels else 0)

# -----------------------------
# Run one Langevin simulation and return step when 4 clusters form
# -----------------------------
def simulate_until_four_clusters(max_steps=200000, check_interval=30, eps=1.7, min_samples=3):
    init_positions = generate_random_particles(n_particles)
    system = ips.LangevinSystem(n_particles, gamma, temperature)

    for i in range(n_particles):
        for d in range(2):
            system.get_positions()[d][i] = init_positions[i][d]
        system.get_velocities()[0][i] = 0.0
        system.get_velocities()[1][i] = 0.0

    sim = ips.IPS_Simulator_Langevin(system)
    pair_force_config = {"type": "LennardJones", "eps": epsilon, "sigma": sigma}
    confinement_config = {"type": "Radial", "rad": rad}
    sim.init(pair_force_config, confinement_config)

    for step in range(max_steps):
        sim.integrate
        if step % check_interval == 0:
            pos = np.array(system.get_positions()).T
            if np.isnan(pos).any():
                print(f"⚠️ NaN detected at step {step}, skipping this run")
                return max_steps  # Treat as failed run
            n_clusters = count_clusters(pos, eps=eps, min_samples=min_samples)
            if n_clusters >= 5:
                return step
    return max_steps

# -----------------------------
# Repeat N simulations and record time to form 4 clusters
# -----------------------------
def repeat_simulations(n_runs=100):
    times = []
    for i in range(n_runs):
        print(f"Run {i+1}/{n_runs}...")
        t = simulate_until_four_clusters()
        times.append(t)
        print(f"  Formed 5 clusters at step {t}")
    return times

# -----------------------------
# Plot histogram
# -----------------------------
def plot_formation_time_distribution(times, label="Simulation"):
    from scipy.stats import poisson
    from sklearn.preprocessing import MinMaxScaler
    from scipy.stats import poisson
    filtered = [t for t in times if t < 200000]
    if not filtered:
        print("No valid data to plot.")
        return

    # Normalize data
    scaler = MinMaxScaler()
    normalized = scaler.fit_transform(np.array(filtered).reshape(-1, 1)).flatten()
    
    # Histogram
    counts, bins, _ = plt.hist(normalized, bins=20, edgecolor='k', alpha=0.6, label='Observed (normalized)')

    # Fit Poisson
    mu = np.mean(filtered)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    poisson_pmf = poisson.pmf(np.round(bin_centers * 20), np.mean(normalized) * 20) * len(filtered), mu) * len(filtered)

    # Overlay Poisson curve
    plt.plot(bin_centers, poisson_pmf, 'r-', lw=2, label=f'Poisson fit (mu={mu:.1f})')

    plt.xlabel("Normalized formation time")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of Cluster Formation Times ({label})")
    plt.legend()
    plt.show()
    print(f"Average formation time (raw): {np.mean(filtered):.2f} steps (out of {len(filtered)} valid runs)")
    print(f"Average formation time: {np.mean(filtered):.2f} steps (out of {len(filtered)} valid runs)")
