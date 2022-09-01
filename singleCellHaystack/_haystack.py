from ._randomization import *
from ._pvalue import *
from ._grid import *
from ._kld import *

from scipy.sparse import isspmatrix

# haystack
# Main function. Can accept AnnData, numpy array and scipy sparse matrices.
# Does not accept numpay matrix objects.
def haystack(x, coord, features=None, scale_coords=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails",
    spline_method="bs", n_randomizations=100, grid_points=None, pseudo=1e-300, verbose=True):

  from anndata import AnnData
  from numpy import ndarray

  if isinstance(x, AnnData) and isinstance(coord, str):
    res = haystack_adata(adata=x, basis=coord, dims=None, scale_coords=scale_coords,
        ngrid_points=ngrid_points, n_genes_to_randomize=n_genes_to_randomize,
        select_genes_randomize_method=select_genes_randomize_method, spline_method=spline_method,
        n_randomizations=n_randomizations, grid_points=grid_points, pseudo=pseudo, verbose=verbose)

  if (isinstance(x, ndarray) or isspmatrix(x)) and isinstance(coord, ndarray):
    res = haystack_array(weights=x, coord=coord, features=features, scale_coords=scale_coords,
        ngrid_points=ngrid_points, n_genes_to_randomize=n_genes_to_randomize,
        select_genes_randomize_method=select_genes_randomize_method, spline_method=spline_method,
        n_randomizations=n_randomizations, grid_points=grid_points, pseudo=pseudo, verbose=verbose)

  return(res)

# haystack_array
# Method for numpy array and scipy sparse matrix objects.
def haystack_array(weights, coord, features=None, scale_coords=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails",
    spline_method="bs", n_randomizations=100, grid_points=None, pseudo=1e-300, verbose=True):

  from pandas import DataFrame
  from sklearn.preprocessing import StandardScaler

  if (verbose):
    print("> entering array method ...")

  exprs = weights

  # Features included?
  if features is None:
    features = np.array(range(exprs.shape[1]))

  # Scale coords.
  if scale_coords:
    if (verbose):
      print("> scaling coordinates ...")

    coord_mean = np.mean(coord, axis=0)
    coord_std = np.std(coord, axis=0)
    coord = (coord - coord_mean) / coord_std

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

  #exprs_sd = np.std(exprs, axis=0)
  sel_zero = exprs_sd == 0
  n_zero = np.sum(sel_zero)
  if (n_zero > 0):
    if (verbose):
      print("> removing", str(n_zero), "genes with zero variance ...")
    sel_nozero = np.invert(sel_zero)
  #   genes = genes[sel_nozero]
    exprs_sd = exprs_sd[sel_nozero] + pseudo
    exprs_mean = exprs_mean[sel_nozero] + pseudo
    exprs = exprs[:, sel_nozero]
    features = features[sel_nozero]

  ncells = exprs.shape[0]
  ngenes = exprs.shape[1]

  # Calculate KLD.
  if grid_points is None:
    grid_points = calculate_grid_points(coord, ngrid_points, verbose=verbose)
  grid_dist = calculate_dist_to_cells(coord, grid_points, verbose=verbose)
  grid_density = calculate_density(grid_dist, verbose=verbose)

  Q = calculate_Q_dist(grid_density, pseudo=pseudo, verbose=verbose)

  KLD = calculate_KLD(grid_density, exprs, Q, verbose=verbose)

  # Calculate CV
  if (verbose):
    print("> calculating feature's CV ...")
  exprs_cv = exprs_sd / exprs_mean

  genes_to_randomize = select_genes_to_randomize(exprs_cv, n_genes_to_randomize, method=select_genes_randomize_method, verbose=verbose)

  # Randomizations.
  KLD_rand = randomize_KLD(grid_density, exprs[:, genes_to_randomize], Q, n_randomizations=n_randomizations, verbose=verbose)

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
    "pval": pval,
    "pval_adj": pval_adj,
    "logpval": logpval,
    "logpval_adj": logpval_adj
  })

  df = df.sort_values("logpval")

  info = {
    "grid_points": grid_points,
    "grid_dist": grid_dist,
    "grid_density": grid_density,
    "Q": Q,
    "KLD_rand": KLD_rand,
    "pval_info": pvalData,
    "genes_to_randomize": genes_to_randomize,
    "exprs_cv": exprs_cv,
  }

  return {
    "results": df,
    "info": info
  }

# haystack_adata
# method for AnnData objects.
def haystack_adata(adata, basis="pca", dims=None, scale_coords=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails", spline_method="bs",
    n_randomizations=100, grid_points=None, pseudo=1e-300, verbose=True):

  if (verbose):
    print("> starting haystack ...")

  if basis in adata.obsm.keys():
    basis_key = basis
  elif f"X_{basis}" in adata.obsm.keys():
    basis_key = f"X_{basis}"

  coord = adata.obsm[basis_key]
  exprs = adata.X
  genes = adata.var_names

  # Check for negative values.
  if (np.sum(exprs < 0).astype(bool)):
    print("_ERROR_ negative values in your data. Hint: pass adata.raw.to_adata() or adata.layers[\"count\"]")
    return(None)

  if dims is not None:
    coord = coord[:, dims]

  res = haystack_array(exprs, coord, features=genes, scale_coords=scale_coords, ngrid_points=ngrid_points,
      n_genes_to_randomize=n_genes_to_randomize, select_genes_randomize_method=select_genes_randomize_method,
      spline_method=spline_method, n_randomizations=n_randomizations, grid_points=grid_points, pseudo=pseudo, verbose=verbose)
  return(res)
