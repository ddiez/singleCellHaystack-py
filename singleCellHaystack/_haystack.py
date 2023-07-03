from ._randomization import *
from ._pvalue import *
from ._grid import *
from ._kld import *
from ._tools import *

from scipy.sparse import isspmatrix

# haystack
# Main function. Can accept AnnData, numpy array and scipy sparse matrices.
# Does not accept numpay matrix objects.
def haystack(x, coord, features=None, layer=None, dims=None, scale_coord=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails", genes_to_randomize=None,
    spline_method="bs", n_randomizations=100, grid_points=None, pseudo=1e-300, random_state=None, verbose=True, kld_method="original"):

  """
  Runs singleCellHaystack.

  :param x: AnnData, numpy array or scipy sparse matrix.
  :param coord: a numpy array with 1D pseudotime, 2D or 3D spatial coordinates or an embedding with any number of dimansions.
  :param features: a list of strings with feature names. If None for AnnData objects is adata.var_names and a numeric index for arrays.
  :param layer: layer to use for AnnData objects. If None then adata.X is used.
  :param scale_coord: whether to scale input coordinates.
  :param ngrid_points: number of grid points.
  :param n_genes_to_randomize: number of genes to use for randomization.
  :param select_genes_randomize_method: method used to select genes for randomization. One of "heavytails" (default) or uniform.
  :param genes_to_randomize: list of genes to randomize.
  :param spline_method: spline method used for randomizations. One of bs (default) or ns.
  :param n_randomizations: number of randomizations.
  :param grid_points: array with grid coordinates.
  :param pseudo: pseudo count added to counts.
  :param random_state: random seed or random state.
  :param verbose: whether to output messages when running haystack.
  :param kld_method: method used to compute KLD.
  :return: A list with singleCellHaystack results.
  :rtype: list

  """
  from anndata import AnnData
  from numpy import ndarray

  res = None

  if isinstance(x, AnnData) and isinstance(coord, str):
    res = haystack_adata(adata=x, basis=coord, layer=layer, dims=dims, scale_coord=scale_coord,
        ngrid_points=ngrid_points, n_genes_to_randomize=n_genes_to_randomize,
        select_genes_randomize_method=select_genes_randomize_method, genes_to_randomize=genes_to_randomize, spline_method=spline_method,
        n_randomizations=n_randomizations, grid_points=grid_points, pseudo=pseudo, random_state=random_state, verbose=verbose, kld_method=kld_method)

  if (isinstance(x, ndarray) or isspmatrix(x)) and isinstance(coord, ndarray):
    res = haystack_array(weights=x, coord=coord, features=features, scale_coord=scale_coord,
        ngrid_points=ngrid_points, n_genes_to_randomize=n_genes_to_randomize,
        select_genes_randomize_method=select_genes_randomize_method, genes_to_randomize=genes_to_randomize, spline_method=spline_method,
        n_randomizations=n_randomizations, grid_points=grid_points, pseudo=pseudo, random_state=random_state, verbose=verbose, kld_method=kld_method)

  if res is None:
    print(
  """
    ERROR: Some of the input data didn't match the expected type:
      * If x is an AnnData object then coord must be a string with the name of the basis to use.
      * If x is a numpy or scipy sparse array then coord must be a numpy array.
  """)
  else:
    return(res)

# haystack_array
# Method for numpy array and scipy sparse matrix objects.
def haystack_array(weights, coord, features=None, scale_coord=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails", genes_to_randomize=None,
    spline_method="bs", n_randomizations=100, grid_points=None, pseudo=1e-300, random_state=None, verbose=True, kld_method="original"):

  from pandas import DataFrame
  from sklearn.preprocessing import StandardScaler
  from numpy.random import RandomState

  if random_state is not None:
    if isinstance(random_state, int):
      random_state = RandomState(random_state)

    if not isinstance(random_state, RandomState):
      print("_ERROR_ invalid random state.")
      return(None)

  if (verbose):
    print("> entering array method ...")

  exprs = weights

  # Features included?
  if features is None:
    features = np.array(range(exprs.shape[1]))

  # Scale coords.
  coord_mean = None
  coord_std = None
  if scale_coord:
    if (verbose):
      print("> scaling coordinates ...")

    coord_mean = np.mean(coord, axis=0)
    coord_std = np.std(coord, axis=0)
    coord = (coord - coord_mean) / coord_std

    #coord, coord_mean, coord_std = scale(coord)

  # Check for negative values.
  if (np.sum(exprs < 0).astype(bool)):
    print("_ERROR_ negative values in your data.")
    return(None)

  # filter genes with zero stdev.
  if (verbose):
    print("> calculating feature stds ...")
  scaler = StandardScaler(with_mean=False)
  scaler_fit = scaler.fit(exprs)
  exprs_sd = np.sqrt(scaler_fit.var_)
  exprs_mean = scaler_fit.mean_

  sel_zero = exprs_sd == 0
  n_zero = np.sum(sel_zero)
  if (n_zero > 0):
    if (verbose):
      print("> removing", str(n_zero), "genes with zero variance ...")
    sel_nozero = np.invert(sel_zero)
    exprs_sd = exprs_sd[sel_nozero] + pseudo
    exprs_mean = exprs_mean[sel_nozero] + pseudo
    exprs = exprs[:, sel_nozero]
    features = features[sel_nozero]

  ncells = exprs.shape[0]
  ngenes = exprs.shape[1]

  # Calculate KLD.
  if grid_points is None:
    grid_points = calculate_grid_points(coord, ngrid_points, random_state=random_state, verbose=verbose)
  grid_dist = calculate_dist_to_cells(coord, grid_points, verbose=verbose)
  grid_density = calculate_density(grid_dist, verbose=verbose)

  Q = calculate_Q_dist(grid_density, pseudo=pseudo, verbose=verbose)
  P = None

  if kld_method == "original":
    KLD = calculate_KLD(grid_density, exprs, Q, verbose=verbose)
  
  if kld_method == "new":
    P = calculate_P_matrix(grid_density, exprs, pseudo=pseudo, verbose=verbose)
    KLD = calculate_KLD2(P, Q, verbose=verbose)

  # Calculate CV
  if (verbose):
    print("> calculating feature's CV ...")
  exprs_cv = exprs_sd / exprs_mean

  if genes_to_randomize is None:
    genes_to_randomize = select_genes_to_randomize(exprs_cv, n_genes_to_randomize, method=select_genes_randomize_method, verbose=verbose)

  # Randomizations.
  if kld_method == "original":
    KLD_rand = randomize_KLD(grid_density, exprs[:, genes_to_randomize], Q, n_randomizations=n_randomizations, random_state=random_state, verbose=verbose)
  
  if kld_method == "new":
    KLD_rand = randomize_KLD2(grid_density, exprs[:, genes_to_randomize], Q, n_randomizations=n_randomizations, random_state=random_state, verbose=verbose)

  # Calculate p.values:
  pvalData = calculate_Pval(KLD, KLD_rand, exprs_cv, exprs_cv[genes_to_randomize], method=spline_method, verbose=verbose)
  pval = pvalData["pval"]
  logpval = pvalData["logpval"]
  #pval_adj = pval * pval.size # Bonferroni correction
  logpval_adj = logpval + np.log10(logpval.size)
  logpval_adj = np.fmin(0, logpval_adj)
  pval_adj = 10 ** logpval_adj

  if (verbose):
    print("> done.")

  # Return results.
  df = DataFrame({
    "gene": features,
    "KLD": KLD,
    "CV": exprs_cv,
    "pval": pval,
    "pval_adj": pval_adj,
    "logpval": logpval,
    "logpval_adj": logpval_adj
  })

  #df = df.sort_values("logpval")

  info = {
    "grid_points": grid_points,
    "grid_dist": grid_dist,
    "grid_density": grid_density,
    "Q": Q,
    "P": P,
    "KLD_rand": KLD_rand,
    "pval_info": pvalData,
    "genes_to_randomize": genes_to_randomize,
    "exprs_cv": exprs_cv,
    "coord_mean": coord_mean,
    "coord_std": coord_std
  }

  res = HaystackResult(result=df, info=info)
  return res

# haystack_adata
# method for AnnData objects.
def haystack_adata(adata, basis="pca", layer=None, dims=None, scale_coord=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails", genes_to_randomize=None, spline_method="bs",
    n_randomizations=100, grid_points=None, pseudo=1e-300, random_state=None, verbose=True, kld_method="original"):

  if (verbose):
    print("> starting haystack ...")

  if basis in adata.obsm.keys():
    basis_key = basis
  elif f"X_{basis}" in adata.obsm.keys():
    basis_key = f"X_{basis}"

  coord = adata.obsm[basis_key]
  if layer is None:
    exprs = adata.X
  else:
    exprs = adata.layers[layer]
  genes = adata.var_names

  # Check for negative values.
  if (np.sum(exprs < 0).astype(bool)):
    print("_ERROR_ negative values in your data. Hint: pass adata.raw.to_adata() or adata.layers[\"count\"]")
    return(None)

  if dims is not None:
    coord = coord[:, :dims]

  res = haystack_array(exprs, coord, features=genes, scale_coord=scale_coord, ngrid_points=ngrid_points,
      n_genes_to_randomize=n_genes_to_randomize, select_genes_randomize_method=select_genes_randomize_method, genes_to_randomize=genes_to_randomize,
      spline_method=spline_method, n_randomizations=n_randomizations, grid_points=grid_points, pseudo=pseudo, random_state=random_state, verbose=verbose, kld_method=kld_method)
  return(res)

class HaystackResult:
  result = None
  info = None

  def __init__(self, result, info):
    self.result=result
    self.info=info

  def top_features(self, n=None, sort_by="logpval_adj"):
    sum = self.result
    sum = sum.sort_values(sort_by)

    if n is not None:
      sum = sum.head(n)
    
    return sum