import numpy as np

from ._backend import is_dask, to_numpy

# x is a matrix of PCA or other embeddings coordinates.
def calculate_grid_points(x, ngrid_points, random_state=None, verbose=False):
  from sklearn.cluster import KMeans

  if (verbose):
    print("> calculating grid points ...")

  # Coordinates are always small (ncells x few dims) compared to the
  # expression matrix, so it is fine to materialize a dask array here:
  # KMeans has no lazy/out-of-core dask support.
  if is_dask(x):
    x = to_numpy(x)

  res = KMeans(n_clusters=ngrid_points, random_state=random_state, n_init=10).fit(x)
  return res.cluster_centers_

def calculate_dist_to_cells(coord, grid_points, verbose=False):
  from scipy.spatial import distance_matrix

  if (verbose):
    print("> calculating distance to cells ...")

  return distance_matrix(to_numpy(coord), to_numpy(grid_points))

# dist is a distance matrix, calculated form grid points and embedding
# coordinates.
def calculate_density(dist, verbose=False):
  if (verbose):
    print("> calculating densities ...")

  bandwidth = np.median(np.min(dist, 1))
  dist_norm = dist / bandwidth
  return np.exp(-dist_norm * dist_norm / 2)
