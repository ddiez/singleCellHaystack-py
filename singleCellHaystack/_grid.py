import numpy as np

# x is a matrix of PCA or other embeddings coordinates.
def calculate_grid_points(x, ngrid_points, random_state=None, verbose=False):
  from sklearn.cluster import KMeans

  if (verbose):
    print("> calculating grid points ...")

  res = KMeans(n_clusters=ngrid_points, random_state=random_state, n_init=10).fit(x)
  return res.cluster_centers_

def calculate_dist_to_cells(coord, grid_points, verbose=False):
  from scipy.spatial import distance_matrix

  if (verbose):
    print("> calculating distance to cells ...")

  return distance_matrix(coord, grid_points)

# dist is a distance matrix, calculated form grid points and embedding
# coordinates.
def calculate_density(dist, verbose=False):
  if (verbose):
    print("> calculating densities ...")

  bandwidth = np.median(np.min(dist, 1))
  dist_norm = dist / bandwidth
  return np.exp(-dist_norm * dist_norm / 2)
