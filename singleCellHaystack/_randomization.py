import numpy as np
from scipy.sparse import isspmatrix
from sklearn.utils import shuffle
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

def randomize_KLD(density, expression, Q, n_randomizations=100, pseudo=1e-300, random_state=None, verbose=False):
  from ._kld import calculate_KLD

  if (verbose):
    print("> calculating randomized KLD ...")
    pbar = tqdm(total=n_randomizations)

  if (isspmatrix(expression)):
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

  if (isspmatrix(expression)):
    expression = expression.tolil()

  ncells = expression.shape[0]
  ngenes = expression.shape[1]

  KLD_rand = np.zeros([ngenes, n_randomizations])

  for n in range(n_randomizations): 
    # TODO: shuffle by column (features) independently.
    expression = shuffle(expression, random_state=random_state)
    #if (isspmatrix(expression)):
    #  for k in range(expression.shape[1]):
    #      expression[:, [k]] = shuffle(expression[:, [k]], random_state=random_state)
    #else:
    #  for k in range(expression.shape[1]):
    #    expression[:, k] = shuffle(expression[:, k], random_state=random_state)
    
    P = calculate_P_matrix(density, expression)
    KLD_rand[:, n] = calculate_KLD2(P, Q)
    if (verbose):
      pbar.update(n=1)

  return KLD_rand
