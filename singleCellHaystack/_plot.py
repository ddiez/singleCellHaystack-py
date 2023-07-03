import matplotlib.pyplot as plt
import numpy as np

def plot_rand_fit(x, type="mean"):
  """
  Plot singleCellHaystack randomization information.

  :param x: singleCellHaystack result.
  :param type: one of "mean" or "sd".

  """

  data = x.info["pval_info"]

  x = data["CV"]
  method = data["method"]

  if type == "mean":
    y = data["rand_mean"]
    y_hat = data["rand_mean_hat"]

    if method == "bs":
      title = "fit: " + method + " n_knots: " + str(data["rand_mean_spline"]["n_knots"]) + " degree: " + str(data["rand_mean_spline"]["degree"])

    if method == "ns":
      title = "fit: " + method + " df: " + str(data["rand_mean_spline"]["df"])

  if type == "sd":
    y = data["rand_sd"]
    y_hat = data["rand_sd_hat"]

    if method == "bs":
      title = "fit: " + method + " n_knots: " + str(data["rand_sd_spline"]["n_knots"]) + " degree: " + str(data["rand_sd_spline"]["degree"])

    if method == "ns":
      title = "fit: " + method + " df: " + str(data["rand_sd_spline"]["df"])

  plt.scatter(x, y, color="black")
  plt.plot(x, y_hat, color="red")
  plt.xlabel("log CV")
  plt.ylabel(type)
  plt.title(title)
  plt.show()

def plot_pval_rank(x):
  """
  Plot singleCellHaystack p.value rank.

  :param x: singleCellHaystack result.

  """

  data = x.info["pval_info"]

  x = np.sort(data["logpval"])

  plt.plot(x, color="black")
  plt.xlabel("Rank")
  plt.ylabel("logpval")
  plt.show()

def plot_pval_hist(x):
  """
  Plot singleCellHaystack p.value histogram.

  :param x: singleCellHaystack result.

  """

  data = x.info["pval_info"]

  x = data["pval"]

  h = plt.hist(x, color="black", bins=30)
  plt.xlabel("P value")
  plt.ylabel("count")
  plt.xlim([0, 1])
  plt.show()

def plot_compare_ranks(res1, res2, sort_by="logpval", xlabel=None, ylabel=None):
  """
  Plot the rank of two singleCellHaystack results.

  This plot facilitates comparing two singleCellHaystack results.

  :param res1: singleCellHaystack result.
  :param res2: singleCellHaystack resuot.
  :param sort_by: what column to use for sorting (rank).
  :param xlabel: label for x-axis (res1).
  :param ylabel: label for y-axis (res2).

  """

  import pandas as pd

  sum1 = res1.result.sort_values(sort_by)
  sum2 = res2.result.sort_values(sort_by)

  if xlabel is None:
    xlabel = sort_by + "_1"

  if ylabel is None:
    ylabel = sort_by + "_2"

  r1 = pd.DataFrame({
    "gene1": sum1.index,
    "rank1": list(range(sum1.index.shape[0])),
  })

  r2 = pd.DataFrame({
    "gene2": sum2.index,
    "rank2": list(range(sum2.index.shape[0])),
  })

  r1 = r1.sort_values("gene1")
  r1 = r1.reset_index(drop=True)

  r2 = r2.sort_values("gene2")
  r2 = r2.reset_index(drop=True)

  d = r1.join(r2)

  plt.plot(d.rank1, d.rank2, "bo", markersize=1)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)

def plot_rand_kld(x, gene=0):
    from scipy.stats import norm
    
    kld_rand = x.info["KLD_rand"]
    X = kld_rand[gene, :]
    X_m = np.mean(X)
    X_s = np.std(X)

    q95 = np.quantile(X, 0.95)
    n95 = norm(loc=X_m, scale=X_s).ppf(0.95)
    
    plt.hist(X, color="black", bins=30)
    plt.grid(False)
    plt.axvline(X_m, color="green", linewidth=1)
    plt.axvline(X_m + X_s, color="purple", linewidth=1)
    plt.axvline(X_m - X_s, color="purple", linewidth=1)
    plt.axvline(q95, color="blue", linewidth=1)
    plt.axvline(n95, color="red", linewidth=1)
    plt.xlabel("KLD_rand")
    plt.ylabel("count")
    plt.show()