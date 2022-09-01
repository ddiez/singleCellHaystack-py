import numpy as np
from scipy.sparse import isspmatrix
from tqdm import tqdm

# density is density as calculated by calculate_density.
def calculate_Q_dist(density, pseudo=1e-300, verbose=False):
  if (verbose):
    print("> calculating Q dist ...")

  Q = np.sum(density, axis=0)
  Q = Q + pseudo
  Q = Q / np.sum(Q)
  return Q

# density (array): ncells x ngrid_points.
# expression (vector): ncells.
# P (vector): ngrid_points.
def calculate_P_dist(density, weights, pseudo=1e-300):
  from anndata import AnnData
  from numpy import ndarray

  if (isspmatrix(weights)):
    index = weights.nonzero()[0]
    P = density[index, :] * weights.data.reshape(-1,1)

  if (isinstance(weights, ndarray)):
    P = density * weights.reshape(-1,1)

  P = np.sum(P, 0)
  P = P + pseudo
  P = P / np.sum(P)
  return P

def calculate_KLD(density, weights, Q, pseudo=1e-300, verbose=False):

  if (isspmatrix(weights)):
    weights = weights.tocsc()

  ngenes = weights.shape[1]

  if (verbose):
    print("> calculating KLD for " + str(ngenes) + " features ...")
    pbar = tqdm(total=ngenes)

  # FIXME: vectorize this computation.
  res = np.zeros(ngenes)
  for k in range(ngenes):
    P = calculate_P_dist(density, weights[:, k])
    res[k] = np.sum(P * np.log(P / Q))
    if (verbose):
      pbar.update(n=1)

  return res
