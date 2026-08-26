import numpy as np

from ._backend import is_sparse, to_numpy

# density is density as calculated by calculate_density.
def calculate_Q_dist(density, pseudo=1e-300, verbose=False):
  if (verbose):
    print("> calculating Q dist ...")

  Q = np.sum(density, axis=0)
  Q = Q + pseudo
  Q = Q / np.sum(Q)
  return Q

# density (array): ncells x ngrid_points.
# weights (vector): ncells.
# P (vector): ngrid_points.
def calculate_P_dist(density, weights, pseudo=1e-300):
  if is_sparse(weights):
    index = weights.nonzero()[0]
    P = density[index, :] * weights.data.reshape(-1, 1)
    P = np.sum(P, axis=0)
  else:
    weights = to_numpy(weights)
    P = np.sum(density * weights.reshape(-1, 1), axis=0)

  P = P + pseudo
  P = P / np.sum(P)
  return P

# density (array): ncells x ngrid_points.
# weights (array): ncells x ngenes.
# P (array): ngenes x ngrid_points.
def calculate_P_matrix(density, weights, pseudo=1e-300, verbose=False):
  if (verbose):
    print("> calculating P for " + str(weights.shape[1]) + " features ...")

  # weights.T @ density sums, for every gene, density rows weighted by that
  # gene's expression -- exactly what calculate_P_dist does per gene, but as
  # a single matrix multiplication. This works unmodified for numpy arrays,
  # scipy sparse matrices and dask arrays (dense or sparse-backed), so no
  # per-backend branching or per-gene python loop is needed here.
  P = to_numpy(weights.T @ density)

  P = P + pseudo
  P = P / np.sum(P, axis=1, keepdims=True)
  return P

def calculate_KLD(density, weights, Q, pseudo=1e-300, verbose=False):
  if (verbose):
    print("> calculating KLD for " + str(weights.shape[1]) + " features ...")

  P = calculate_P_matrix(density, weights, pseudo=pseudo, verbose=False)
  return calculate_KLD2(P, Q, verbose=verbose)

def calculate_KLD2(P, Q, verbose=False):
  if (verbose):
    print("> calculating KLD ...")

  kld = np.sum(P * np.log(P / Q), axis=1)
  return kld
