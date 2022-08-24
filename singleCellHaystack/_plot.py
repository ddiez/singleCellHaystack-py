import matplotlib.pyplot as plt
import numpy as np

def plot_rand_fit(x, type="mean"):
  data = x["info"]["pval_info"]

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
  data = x["info"]["pval_info"]

  x = np.sort(data["logpval"])

  plt.plot(x, color="black")
  plt.xlabel("Rank")
  plt.ylabel("logpval")
  plt.show()

def plot_pval_hist(x):
  data = x["info"]["pval_info"]

  x = data["pval"]

  h = plt.hist(x, color="black", bins=30)
  plt.xlim([0, 1])
  plt.show()

def plot_compare_ranks(res1, res2, sort_by="logpval", xlabel=None, ylabel=None):
  import pandas as pd

  sum1 = res1["results"].sort_values(sort_by)
  sum2 = res2["results"].sort_values(sort_by)

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
