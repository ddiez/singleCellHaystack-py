import numpy as np
from numpy import ndarray

import pandas as pd
from anndata import AnnData

from scipy.spatial import distance_matrix
from random import sample
import scipy.sparse as sp
#from sklearn.preprocessing import StandardScaler
from scipy.sparse import isspmatrix

from tqdm import tqdm

# def haystack(x, basis="pca", method="highD", dims=None, ngrid_points=100,
#     pseudo=1e-300, n_genes_to_randomize=100, n_randomizations=100, grid_points=None):
#
#  from numpy import ndarray
#
#
#   if isinstance(x, AnnData) and isinstance(basis, str):
#     haystack_adata(x, method=method, basis=basis, dims=None, ngrid_points=100,
#         pseudo=1e-300, n_genes_to_randomize=100, n_randomizations=100, grid_points=None)
#   if isinstance(x, ndarray) and isinstance(basis, ndarray):
#     haystack_ndarray(weights=x, coord=basis)
#
# def haystack_ndarray(weights: ndarray, coord: ndarray):
#   print("ndarray method")
#   return(0)
#
def haystack_sparse(exprs: ndarray, coord: ndarray, scale_coords=True, ngrid_points=100,
    n_genes_to_randomize=100, select_genes_randomize_method="heavytails",
    spline_method="bs", n_randomizations=100, grid_points=None, verbose=True, pseudo=1e-300):

  if (verbose):
    print("> entering sparse method ...")

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
  from sklearn.preprocessing import StandardScaler
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

  ncells = exprs.shape[0]
  ngenes = exprs.shape[1]

  # Calculate KLD.
  if grid_points is None:
    grid_points = calculate_grid_points(coord, ngrid_points, verbose=verbose)
  grid_dist = calculate_dist_to_cells(coord, grid_points, verbose=verbose)
  grid_density = calculate_density(grid_dist, verbose=verbose)

  Q = calculate_Q_dist(grid_density, pseudo=pseudo, verbose=verbose)

  KLD = calculate_KLD_sparse(grid_density, exprs, Q, verbose=verbose)

  # Calculate CV
  if (verbose):
    print("> calculating feature means ...")
  #exprs_mean = scaler_fit.mean_[sel_nozero]
  #exprs_mean = np.mean(exprs, axis=0) + pseudo
  exprs_cv = exprs_sd / exprs_mean

  genes_to_randomize = select_genes_to_randomize(exprs_cv, n_genes_to_randomize, method=select_genes_randomize_method, verbose=verbose)

  # Randomizations.
  KLD_rand = randomize_KLD_sparse(grid_density, exprs[:, genes_to_randomize], Q, n_randomizations=n_randomizations, verbose=verbose)

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
  df = pd.DataFrame({
    #"gene": genes,
    "KLD": KLD,
    "pval": pval,
    "pval_adj": pval_adj,
    "logpval": logpval,
    "logpval_adj": logpval_adj
  })

  df = df.sort_values("logpval_adj")

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

# def haystack_adata(adata, method="highD", basis="pca", dims=None, ngrid_points=100,
#     pseudo=1e-300, n_genes_to_randomize=100, n_randomizations=100, grid_points=None):
#  print("anndata method")
#  return(0)

def haystack(adata, basis="pca", dims=None, scale_coords=True, ngrid_points=100,
    pseudo=1e-300, n_genes_to_randomize=100, select_genes_randomize_method="heavytails",
    spline_method="bs", n_randomizations=100, grid_points=None, verbose=True):

  if (verbose):
    print("> starting haystack ...")

  if basis in adata.obsm.keys():
    basis_key = basis
  elif f"X_{basis}" in adata.obsm.keys():
    basis_key = f"X_{basis}"

  coord = adata.obsm[basis_key]
  exprs = adata.X
  genes = adata.var_names

  # Scale coords.
  if scale_coords:
    if (verbose):
      print("> scaling coordinates ...")

    coord_mean = np.mean(coord, axis=0)
    coord_std = np.std(coord, axis=0)
    coord = (coord - coord_mean) / coord_std

  # FIXME: workaround, fix properly.
  if isspmatrix(exprs):
    if (verbose):
      print("> converting to dense array ...")

    exprs = exprs.toarray()

  # Check for negative values.
  if (np.sum(exprs < 0).astype(bool)):
    print("_ERROR_ negative values in your data. Hint: pass adata.raw.to_adata() or adata.layers[\"count\"]")
    return(None)

  if dims is not None:
    coord = coord[:, dims]

  # We could call a haystack matrix based function here.
  # haystack_array(coord, exprs, genes, ...)

  # filter genes with zero stdev.
  if (verbose):
    print("> calculating feature stds ...")
  exprs_sd = np.std(exprs, axis=0)
  sel_zero = exprs_sd == 0
  n_zero = np.sum(sel_zero)
  if (n_zero > 0):
    if (verbose):
      print("> removing", str(n_zero), "genes with zero variance ...")
    sel_nozero = np.invert(sel_zero)
    genes = genes[sel_nozero]
    exprs_sd = exprs_sd[sel_nozero]
    exprs = exprs[:, sel_nozero]

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
    print("> calculating feature means ...")
  exprs_mean = np.mean(exprs, axis=0) + pseudo
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
  df = pd.DataFrame({
    "gene": genes,
    "KLD": KLD,
    "pval": pval,
    "pval_adj": pval_adj,
    "logpval": logpval,
    "logpval_adj": logpval_adj
  })

  df = df.sort_values("logpval_adj")

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

# x is a matrix of PCA or other embeddings coordinates.
def calculate_grid_points(x, ngrid_points, random_state=None, verbose=False):
  from sklearn.cluster import KMeans

  if (verbose):
    print("> calculating grid points ...")

  res = KMeans(n_clusters=ngrid_points, random_state=random_state).fit(x)
  return res.cluster_centers_

def calculate_dist_to_cells(coord, grid_points, verbose=False):
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
def calculate_P_dist(density, expression, pseudo=1e-300):
  P = density * expression.reshape(-1,1)

  P = np.sum(P, 0)
  P = P + pseudo
  P = P / np.sum(P)
  return P

def calculate_P_dist_sparse(density, expression, pseudo=1e-300):

  index = expression.nonzero()[0]
  P = density[index, :] * expression.data.reshape(-1,1)

  P = np.sum(P, 0)
  P = P + pseudo
  P = P / np.sum(P)
  return P

def calculate_KLD(density, expression, Q, pseudo=1e-300, verbose=False):

  ngenes = expression.shape[1]

  if (verbose):
    print("> calculating KLD for " + str(ngenes) + " genes ...")
    pbar = tqdm(total=ngenes)

  # FIXME: vectorize this computation.
  res = np.zeros(ngenes)
  for k in range(ngenes):
    P = calculate_P_dist(density, expression[:, k])
    res[k] = np.sum(P * np.log(P / Q))
    if (verbose):
      pbar.update(n=1)

  return res

def calculate_KLD_sparse(density, expression, Q, pseudo=1e-300, verbose=False):
  expression = expression.tocsc()
  ngenes = expression.shape[1]

  if (verbose):
    print("> calculating KLD for " + str(ngenes) + " genes ...")
    pbar = tqdm(total=ngenes)

  # FIXME: vectorize this computation.
  res = np.zeros(ngenes)
  for k in range(ngenes):
    P = calculate_P_dist_sparse(density, expression[:, [k]])
    res[k] = np.sum(P * np.log(P / Q))
    if (verbose):
      pbar.update(n=1)

  return res


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
  if (verbose):
    print("> calculating randomized KLD ...")
    pbar = tqdm(total=n_randomizations)

  ncells=expression.shape[0]
  ngenes=expression.shape[1]

  KLD_rand=np.zeros([ngenes, n_randomizations])

  for n in range(n_randomizations):
    shuffled_cells = sample(range(ncells), ncells)
    KLD_rand[:, n] = calculate_KLD(density, expression[shuffled_cells, :], Q)
    if (verbose):
      pbar.update(n=1)

  return KLD_rand

def randomize_KLD_sparse(density, expression, Q, n_randomizations=100, pseudo=1e-300, verbose=False):
  if (verbose):
    print("> calculating randomized KLD ...")
    pbar = tqdm(total=n_randomizations)

  ncells=expression.shape[0]
  ngenes=expression.shape[1]

  KLD_rand=np.zeros([ngenes, n_randomizations])

  for n in range(n_randomizations):
    shuffled_cells = sample(range(ncells), ncells)
    KLD_rand[:, n] = calculate_KLD_sparse(density, expression[shuffled_cells, :], Q)
    if (verbose):
      pbar.update(n=1)

  return KLD_rand

def estimate_spline_param(x, y, method="bs"):
  from sklearn.preprocessing import SplineTransformer
  from sklearn.preprocessing import FunctionTransformer
  from sklearn.linear_model import LinearRegression
  from sklearn.pipeline import Pipeline
  from sklearn.model_selection import GridSearchCV
  from patsy import cr

  if method == "bs":
    pip = Pipeline([
      ["transformer", SplineTransformer()],
      ["estimator", LinearRegression()]
    ])

    min_knots = 2
    max_knots = 10
    min_degree = 1
    max_degree = 5
    param = {
      "transformer__n_knots": list(range(min_knots, max_knots)),
      "transformer__degree": list(range(min_degree, max_degree))
    }

    cv = GridSearchCV(pip, param, cv=10)
    cv_res = cv.fit(x, y)

    info = {
      "method": "bs",
      "n_knots": cv_res.best_params_["transformer__n_knots"],
      "degree": cv_res.best_params_["transformer__degree"]
    }

  if method == "ns":
    NaturalSplineTransformer = FunctionTransformer(cr)
    pip = Pipeline([
      ["transformer", NaturalSplineTransformer],
      ["estimator", LinearRegression()]
    ])

    # Set dict of degrees of freedom for CV.
    min_df = 3
    max_df = 10
    df = []
    for k in range(min_df, max_df):
      df.append({"df": k})

    param = {
      "transformer__kw_args": df
    }

    cv = GridSearchCV(pip, param, cv=10)
    cv_res = cv.fit(x, y)

    info = {
      "method": "ns",
      "df": cv_res.best_params_["transformer__kw_args"]["df"]
    }

  return info

def calculate_KLD_fit(x, y, method="bs"):
  from sklearn.preprocessing import SplineTransformer
  from sklearn.preprocessing import FunctionTransformer
  from sklearn.linear_model import LinearRegression
  from sklearn.pipeline import Pipeline
  from patsy import cr

  info = estimate_spline_param(x, y, method=method)
  if method == "bs":
    n_knots = info["n_knots"]
    degree = info["degree"]

    pip = Pipeline([
      ["transformer", SplineTransformer(n_knots=n_knots, degree=degree)],
      ["estimator", LinearRegression()]
    ])
    model = pip.fit(x, y)
    y_hat = model.predict(x)

  if method == "ns":
    df = info["df"]

    NaturalSplineTransformer = FunctionTransformer(cr, kw_args={"df": df})
    pip = Pipeline([
      ["transformer", NaturalSplineTransformer],
      ["estimator", LinearRegression()]
    ])
    model = pip.fit(x, y)
    y_hat = model.predict(x)

  return {
    "model": model,
    "spline": info,
    "y_hat": y_hat
  }

def calculate_Pval(KLD, KLD_rand, cv, cv_rand, method="bs", verbose=False):
  from sklearn.preprocessing import SplineTransformer
  from sklearn.linear_model import LinearRegression
  from sklearn.pipeline import Pipeline
  from patsy import cr
  from scipy.stats import norm

  if (verbose):
    print("> calculating P values ...")

  KLD_log = np.log(KLD)
  KLD_rand_log = np.log(KLD_rand)
  KLD_rand_mean = np.mean(KLD_rand_log, axis=1)
  KLD_rand_sd = np.std(KLD_rand_log, axis=1)
  cv_log = np.log(cv).reshape(-1, 1)
  cv_rand_log = np.log(cv_rand).reshape(-1, 1)

  KLD_rand_mean_fit = calculate_KLD_fit(cv_rand_log, KLD_rand_mean, method=method)
  KLD_rand_sd_fit = calculate_KLD_fit(cv_rand_log, KLD_rand_sd, method=method)

  KLD_rand_mean_model = KLD_rand_mean_fit["model"]
  KLD_rand_sd_model = KLD_rand_sd_fit["model"]

  KLD_mean = KLD_rand_mean_model.predict(cv_log)
  KLD_sd = KLD_rand_sd_model.predict(cv_log)

  logpval = norm.logsf(KLD_log, loc=KLD_mean, scale=KLD_sd)/np.log(10)
  pval = 10 ** logpval

  return {
    "pval": pval,
    "logpval": logpval,
    "method": method,
    "CV": cv_rand_log,
    "rand_mean": KLD_rand_mean,
    "rand_sd": KLD_rand_sd,
    "rand_mean_model": KLD_rand_mean_model,
    "rand_sd_model": KLD_rand_sd_model,
    "rand_mean_spline": KLD_rand_mean_fit["spline"],
    "rand_sd_spline": KLD_rand_sd_fit["spline"],
    "rand_mean_hat": KLD_rand_mean_fit["y_hat"],
    "rand_sd_hat": KLD_rand_sd_fit["y_hat"]
  }
