import numpy as np
from scipy.sparse import issparse

# Small helpers that let the rest of the package treat numpy arrays, scipy
# sparse matrices and dask arrays uniformly, instead of scattering
# isinstance/isspmatrix checks (and now dask checks) across every function.
# dask is an optional dependency: is_dask() returns False when it isn't
# installed rather than raising.

def is_sparse(x):
  """True if x is a scipy sparse matrix/array."""
  return issparse(x)

def is_pydata_sparse(x):
  """True if x is a pydata/sparse array (the format dask uses for sparse-backed arrays)."""
  try:
    from sparse import SparseArray
  except ImportError:
    return False
  return isinstance(x, SparseArray)

def is_dask(x):
  """True if x is a dask array."""
  try:
    from dask.array import Array
  except ImportError:
    return False
  return isinstance(x, Array)

def to_numpy(x):
  """Materialize x (dask array, scipy/pydata sparse array, or array-like) into a plain numpy array."""
  if is_dask(x):
    x = x.compute()
  if is_sparse(x):
    x = x.toarray()
  if is_pydata_sparse(x):
    x = x.todense()
  return np.asarray(x)

def calculate_mean_std(x, axis=0, verbose=False):
  """
  Calculate mean and (population) standard deviation of x along axis, without
  forcing a specific backend (numpy, scipy sparse or dask are all supported).
  """
  if (verbose):
    print("> calculating mean and standard deviation ...")

  if is_sparse(x):
    sq = x.multiply(x)
  else:
    sq = x ** 2

  mean = to_numpy(x.mean(axis=axis)).ravel()
  mean_sq = to_numpy(sq.mean(axis=axis)).ravel()

  # clip to avoid tiny negative variances caused by floating point error.
  var = np.clip(mean_sq - mean ** 2, a_min=0, a_max=None)
  std = np.sqrt(var)

  return mean, std
