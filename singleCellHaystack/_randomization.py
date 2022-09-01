import numpy as np
from scipy.sparse import isspmatrix
from random import sample
from tqdm import tqdm

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

def randomize_KLD(density, expression, Q, n_randomizations=100, pseudo=1e-300, verbose=False):
  from ._haystack import calculate_KLD

  if (verbose):
    print("> calculating randomized KLD ...")
    pbar = tqdm(total=n_randomizations)

  if (isspmatrix(expression)):
    expression = expression.tocsc()

  ncells=expression.shape[0]
  ngenes=expression.shape[1]

  KLD_rand=np.zeros([ngenes, n_randomizations])

  for n in range(n_randomizations):
    shuffled_cells = sample(range(ncells), ncells)
    KLD_rand[:, n] = calculate_KLD(density, expression[shuffled_cells, :], Q)
    if (verbose):
      pbar.update(n=1)

  return KLD_rand
