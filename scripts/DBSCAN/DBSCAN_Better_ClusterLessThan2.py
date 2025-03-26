import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

positions_final = np.array([p.get_positions()[0], p.get_positions()[1]]).T

# Compute eps: based on k-distance curve + maximum curvature
def find_eps_auto(positions, k=3):
    neigh = NearestNeighbors(n_neighbors=k)
    distances, _ = neigh.fit(positions).kneighbors(positions)
    k_distances = distances[:, -1]
    sorted_distances = np.sort(k_distances)[::-1]
    grad = np.gradient(sorted_distances)
    curvature = np.gradient(grad)
    idx = np.argmax(curvature)
    eps = sorted_distances[idx]
    return sorted_distances, idx, eps

# Set parameters based on dimensionality
dimension = positions_final.shape[1]
k = 2 * dimension - 1
min_samples = k + 1

# Compute eps and perform DBSCAN
sorted_distances, idx, eps = find_eps_auto(positions_final, k=k)
db = DBSCAN(eps=eps, min_samples=min_samples)
labels = db.fit_predict(positions_final)
sil_score = silhouette_score(positions_final, labels) if len(np.unique(labels)) > 1 else 0

plt.figure(figsize=(15, 5))

# K-Distance curve plot
plt.subplot(1, 2, 1)
plt.plot(sorted_distances, label='k-distance curve', color='orange')
plt.axvline(x=idx, color='red', linestyle='--', label=f'eps = {eps:.2f}')
plt.xlabel("Sorted Points (descending)")
plt.ylabel(f"{k}th nearest neighbor distance")
plt.title("K-Distance Graph with Auto-selected Eps")
plt.legend()

# DBSCAN clustering visualization
plt.subplot(1, 2, 2)
for label in np.unique(labels):
    mask = labels == label
    color = 'gray' if label == -1 else plt.cm.tab10(label % 10)
    name = "Noise" if label == -1 else f"Cluster {label}"
    plt.scatter(positions_final[mask, 0], positions_final[mask, 1], c=[color], edgecolor='k', label=name)
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'DBSCAN Clustering (eps={eps:.2f}, min_samples={min_samples})\nSilhouette: {sil_score:.2f}')
plt.legend()
plt.tight_layout()
plt.show()
