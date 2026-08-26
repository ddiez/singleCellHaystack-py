import numpy as np
from sklearn.utils import shuffle, check_random_state
from tqdm import tqdm

from ._backend import is_sparse

def select_genes_to_randomize(x, ngenes=100, method="heavytails", tail=10, verbose=False):
  if (verbose):
    print("> selecting genes to randomize ...")

  index = np.argsort(x)

  if (method == "uniform"):
    return index[np.linspace(0, index.shape[0]-1, ngenes, dtype="int")]

  if (method == "heavytails"):
    ls = index[:tail]
    rs = index[-tail:]
    ngenes = ngenes - tail * 2
    ms = index[np.linspace(tail, index.shape[0]-tail-1, ngenes, dtype="int")]
    return np.concatenate([ls, ms, rs])

def randomize_KLD(density, expression, Q, n_randomizations=100, pseudo=1e-300, random_state=None, verbose=False):
  from ._kld import calculate_KLD

  if (verbose):
    print("> calculating randomized KLD ...")
    pbar = tqdm(total=n_randomizations)

  # lil is faster than csr/csc for repeated row slicing below.
  if is_sparse(expression):
    expression = expression.tolil()

  ncells = expression.shape[0]
  ngenes = expression.shape[1]

  KLD_rand = np.zeros([ngenes, n_randomizations])

  for n in range(n_randomizations):
    shuffled_cells = shuffle(range(ncells), random_state=random_state)
    KLD_rand[:, n] = calculate_KLD(density, expression[shuffled_cells, :], Q)
    if (verbose):
      pbar.update(n=1)

  return KLD_rand

def randomize_KLD2(density, expression, Q, n_randomizations=100, pseudo=1e-300, random_state=None, verbose=False):
  from ._kld import calculate_KLD2
  from ._kld import calculate_P_matrix

  if (verbose):
    print("> calculating randomized KLD ...")
    pbar = tqdm(total=n_randomizations)

  ncells = expression.shape[0]
  ngenes = expression.shape[1]
  rng = check_random_state(random_state)

  KLD_rand = np.zeros([ngenes, n_randomizations])

  for n in range(n_randomizations):
    # Draw a fresh permutation of the (original) cells for every
    # randomization and index with it directly: this works the same way for
    # numpy arrays, scipy sparse matrices and dask arrays, so it doesn't need
    # any backend-specific handling.
    shuffled_cells = rng.permutation(ncells)
    P = calculate_P_matrix(density, expression[shuffled_cells, :], pseudo=pseudo)
    KLD_rand[:, n] = calculate_KLD2(P, Q)
    if (verbose):
      pbar.update(n=1)

  return KLD_rand
