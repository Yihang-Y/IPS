import numpy as np
from sklearn.neighbors import NearestNeighbors
from kneed import KneeLocator
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from matplotlib.cm import tab20

# Retrieve final particle positions
positions_final = np.array([p.get_positions()[0], p.get_positions()[1]]).T

# Core Algorithm
def calculate_smoothness(sorted_distances):
    """Evaluate smoothness of the last 30% of the curve"""
    tail_length = int(len(sorted_distances) * 0.3)
    tail = sorted_distances[-tail_length:]
    grad = np.gradient(np.gradient(tail))
    return np.std(grad)

def auto_interp_method(sorted_distances, smooth_threshold=0.05):
    """Automatically select interpolation method"""
    smoothness = calculate_smoothness(sorted_distances)
    return 'polynomial' if smoothness > smooth_threshold else 'interp1d'

def enhanced_kneedle(positions, k=5):
    """Enhanced elbow detection algorithm"""
    # Calculate k-distance curve
    neigh = NearestNeighbors(n_neighbors=k)
    nbrs = neigh.fit(positions)
    distances, _ = nbrs.kneighbors(positions)
    k_distances = np.sort(distances[:, -1])[::-1]
    
    # Adaptive parameter configuration
    interp_method = auto_interp_method(k_distances)
    method_params = {'interp_method': interp_method}
    if interp_method == 'polynomial':
        method_params['polynomial_degree'] = 3
    
    # Perform elbow detection
    x = np.arange(len(k_distances))
    kneedle = KneeLocator(
        x, k_distances,
        curve='convex', 
        direction='decreasing',
        online=True,
        **method_params
    )
    
    # Exception handling strategy
    if kneedle.elbow is None or kneedle.elbow < len(k_distances)*0.1:
        return enhanced_kneedle(positions, k+1)  # Retry with incremented k
        
    return k_distances, kneedle.elbow, k_distances[kneedle.elbow]

# Automatic Parameter Generation
def auto_params(dimension):
    """Dimension-aware parameter generation (2D → k=3, min_samples=4)"""
    k = 2 * dimension - 1
    return k, k + 1

# Main Pipeline
dimension = positions_final.shape[1]  # Dimension detection (2D)
k, min_samples = auto_params(dimension)

# Execute algorithm
try:
    sorted_distances, idx, eps = enhanced_kneedle(positions_final, k)
except ValueError as e:
    print(f"Error: {e}")
    exit()

# Visualization System
plt.figure(figsize=(18, 6))

# 1: Adaptive K-Distance Analysis
plt.subplot(1, 3, 1)
plt.plot(sorted_distances, color='navy', linewidth=2, label='k-distance curve')
plt.axvline(idx, color='crimson', linestyle='--', 
            label=f'Adaptive eps={eps:.2f} (index {idx})')
plt.axvspan(0, idx, alpha=0.1, color='red', label='Noise Region')
plt.axvspan(idx, len(sorted_distances), alpha=0.1, color='green', label='Cluster Region')
plt.xlabel("Sorted Points (descending)")
plt.ylabel(f"{k}th NN Distance")
plt.title(f"Enhanced K-Distance Analysis\n(Dimension={dimension})")
plt.grid(True, alpha=0.3)
plt.legend()

# 2: Clustering Results
plt.subplot(1, 3, 2)
db = DBSCAN(eps=eps, min_samples=min_samples)
labels = db.fit_predict(positions_final)

# Calculate quality metrics
sil_score = silhouette_score(positions_final, labels) if len(np.unique(labels))>1 else 0

# Modified color generation section
for label in np.unique(labels):
    mask = labels == label
    if label == -1:
        color = [0.7, 0.7, 0.7]  # Noise remains gray
        marker = 'o'
    else:
        color = plt.cm.tab20(label % 20)  # Use tab20 colormap
        marker = ['o', '^', 's', 'D'][label % 4]  # Cycle through 4 markers
    
    plt.scatter(positions_final[mask,0], positions_final[mask,1], 
                c=[color], 
                s=50, 
                edgecolor='k', 
                alpha=0.8,
                marker=marker,
                label=f'Cluster {label}' if label!=-1 else 'Noise')

plt.title(f'Clustering Result\neps={eps:.2f}, min_samples={min_samples}\nSilhouette: {sil_score:.2f}')
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.grid(True, alpha=0.3)

# 3: Smoothness Analysis
plt.subplot(1, 3, 3)
window_size = 50
smoothness = [calculate_smoothness(sorted_distances[i:i+window_size]) 
              for i in range(len(sorted_distances)-window_size)]
plt.plot(smoothness, color='purple', label='Local Smoothness')
plt.axhline(0.05, color='orange', linestyle='--', label='Threshold')
plt.xlabel('Sliding Window Position')
plt.ylabel('Smoothness Level')
plt.title('Local Curve Smoothness Analysis')
plt.legend()

plt.tight_layout()
plt.show()

print(f"[Success] Parameters: k={k}, eps={eps:.2f}, min_samples={min_samples}")
print(f"          Noise Points: {np.sum(labels==-1)}/{len(labels)}")
print(f"          Cluster Quality: Silhouette = {sil_score:.2f}")